# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Build a pre-equilibrated periodic **solvent box** for a single CHARMM residue, for use
as a ``kind: box`` entry in the PDB repository's ``solvent`` collection.  Such a box is
what VMD's ``solvate`` plugin needs to bulk-solvate with anything other than its built-in
TIP3P water (``-spsf``/``-spdb``/``-ws``/``-ks``); see
``docs/design/solvent-collection.md``.

This module holds the deterministic, dependency-light core -- box geometry and cubic
packing -- separately from the psfgen/NAMD build orchestration, so the geometry can be
unit-tested without a force field or MD engine.
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)

# Avogadro's number (1/mol); 1 cm^3 = 1e24 A^3
_N_AVOGADRO = 6.02214076e23
_CM3_PER_A3 = 1.0e-24


def box_edge_for_density(nmol: int, molweight: float, density_gcc: float) -> float:
    """
    Edge length (Å) of a cube holding ``nmol`` molecules of molar mass ``molweight``
    (g/mol) at bulk density ``density_gcc`` (g/cm^3).

    Volume = (nmol * molweight / N_A) / density, converted from cm^3 to Å^3.  This sizes
    the *initial* box; the NPT equilibration relaxes it to the true equilibrium edge,
    which is what gets recorded for ``solvate -ws``.
    """
    if nmol <= 0 or molweight <= 0 or density_gcc <= 0:
        raise ValueError('nmol, molweight, and density_gcc must all be positive')
    mass_g = nmol * molweight / _N_AVOGADRO
    volume_cm3 = mass_g / density_gcc
    volume_a3 = volume_cm3 / _CM3_PER_A3
    return volume_a3 ** (1.0 / 3.0)


def _random_rotations(n: int, rng: np.random.Generator) -> np.ndarray:
    """Return ``n`` uniformly-random 3x3 rotation matrices (via random unit quaternions)."""
    # Marsaglia/Shoemake: uniform quaternions -> rotation matrices, no scipy dependency
    u1, u2, u3 = rng.random(n), rng.random(n), rng.random(n)
    q = np.empty((n, 4))
    q[:, 0] = np.sqrt(1 - u1) * np.sin(2 * np.pi * u2)
    q[:, 1] = np.sqrt(1 - u1) * np.cos(2 * np.pi * u2)
    q[:, 2] = np.sqrt(u1) * np.sin(2 * np.pi * u3)
    q[:, 3] = np.sqrt(u1) * np.cos(2 * np.pi * u3)
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    R = np.empty((n, 3, 3))
    R[:, 0, 0] = 1 - 2 * (y * y + z * z); R[:, 0, 1] = 2 * (x * y - z * w); R[:, 0, 2] = 2 * (x * z + y * w)
    R[:, 1, 0] = 2 * (x * y + z * w); R[:, 1, 1] = 1 - 2 * (x * x + z * z); R[:, 1, 2] = 2 * (y * z - x * w)
    R[:, 2, 0] = 2 * (x * z - y * w); R[:, 2, 1] = 2 * (y * z + x * w); R[:, 2, 2] = 1 - 2 * (x * x + y * y)
    return R


def cubic_lattice_sites(nmol: int, edge: float) -> np.ndarray:
    """
    ``nmol`` cell-centered sites of the smallest cubic lattice that holds them, spanning
    an ``edge``-length cube.  Returns an ``(nmol, 3)`` array of centers in ``[0, edge)``.
    """
    ncell = int(np.ceil(round(nmol ** (1.0 / 3.0), 6)))
    while ncell ** 3 < nmol:
        ncell += 1
    spacing = edge / ncell
    grid = (np.arange(ncell) + 0.5) * spacing
    sites = np.array([(x, y, z) for x in grid for y in grid for z in grid])
    return sites[:nmol]


def pack_cubic(mol_coords: np.ndarray, nmol: int, edge: float, seed=None,
               jitter: float = 0.0) -> np.ndarray:
    """
    Place ``nmol`` randomly-oriented copies of a molecule (whose atom coordinates are
    ``mol_coords``, an ``(natom, 3)`` array) on the cubic lattice of
    :func:`cubic_lattice_sites`.

    Each copy is centered on its atoms' centroid, given a uniformly-random orientation,
    and translated to a lattice site (plus optional ``jitter``).  Returns an
    ``(nmol, natom, 3)`` array of placed coordinates.  Overlaps left for the downstream
    NPT relaxation to resolve.
    """
    rng = np.random.default_rng(seed)
    centered = mol_coords - mol_coords.mean(axis=0)
    sites = cubic_lattice_sites(nmol, edge)
    rots = _random_rotations(nmol, rng)
    placed = np.einsum('nij,aj->nai', rots, centered)   # rotate each copy
    if jitter:
        sites = sites + rng.uniform(-jitter, jitter, size=sites.shape)
    placed += sites[:, None, :]
    return placed


def _segid(i: int) -> str:
    """A 4-character base-36 segid for the i-th molecule (covers 36^4 = 1.6M copies)."""
    digits = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    s = ''
    for _ in range(4):
        i, r = divmod(i, 36)
        s = digits[r] + s
    return s


def write_box_pdb(placed: np.ndarray, template_lines: list, output_pdb: str) -> list:
    """
    Write an ``nmol``-molecule box PDB by rewriting one molecule's ATOM lines
    (``template_lines``, ``natom`` fixed-format PDB records) once per placed copy, with a
    fresh serial, a per-molecule ``resSeq`` of 1, the placed coordinates, and a unique
    4-character ``segID`` per copy (so psfgen can build one segment per molecule).

    ``placed`` is the ``(nmol, natom, 3)`` array from :func:`pack_cubic`.  Returns the
    list of segids used, in order.
    """
    nmol, natom, _ = placed.shape
    if len(template_lines) != natom:
        raise ValueError(f'template has {len(template_lines)} atom lines but placed has {natom} atoms')
    tmpl = [ln.rstrip('\n').ljust(80) for ln in template_lines]
    segids = []
    serial = 0
    with open(output_pdb, 'w') as f:
        for i in range(nmol):
            sid = _segid(i)
            segids.append(sid)
            for a in range(natom):
                serial += 1
                x, y, z = placed[i, a]
                ln = tmpl[a]
                rec = (ln[:6]
                       + f'{serial % 100000:5d}'
                       + ln[11:22]
                       + f'{1:4d}'                      # resSeq = 1 for this molecule
                       + ln[26:30]
                       + f'{x:8.3f}{y:8.3f}{z:8.3f}'
                       + ln[54:72]
                       + f'{sid:<4s}'
                       + ln[76:80])
                f.write(rec.rstrip() + '\n')
        f.write('END\n')
    return segids


def atom_lines_of(pdb_path: str) -> list:
    """Return the ATOM/HETATM record lines of a PDB file (in order)."""
    with open(pdb_path) as f:
        return [ln for ln in f if ln.startswith(('ATOM', 'HETATM'))]


def coords_of(atom_lines: list) -> np.ndarray:
    """Parse the (natom, 3) coordinate array from fixed-format PDB ATOM lines."""
    return np.array([[float(ln[30:38]), float(ln[38:46]), float(ln[46:54])] for ln in atom_lines])
