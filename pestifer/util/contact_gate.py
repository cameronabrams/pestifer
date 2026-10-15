# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Refuse to hand NAMD a system with an unresolvable steric overlap.

A RATTLE constraint failure ends a build outright, and NAMD reports it as an atom serial
after the run has started -- typically hours in, on a cluster, with the diagnosis left to
whoever reads the log.  The overlap that caused it was present in the coordinates before
NAMD ever launched and is cheap to find there.

This is the same bargain ``NAMDScripter._verify_params_cover_psf`` already makes for force-field
parameters, for the same reason: a missing parameter is fatal to NAMD regardless, so the only
question is whether the failure happens somewhere it can be diagnosed.  A steric overlap is not
*always* fatal -- a minimization is exactly the thing that fixes one -- so this gates **dynamics
only**, never a minimize.  Gating a minimize would make an overlapped system unfixable.

**What counts as an overlap.**  Heavy atoms only: X-H distances are what RATTLE constrains, and
a legitimate H...H contact reaches ~1.8 A, so including hydrogens would mean choosing between
false positives and a useless threshold.  Pairs that are 1-2 or 1-3 bonded are excluded, which
is what makes this work across a polymer: a peptide C-N bond is 1.33 A and the CA...N across it
is ~2.4 A, and neither is a clash.  Anything else below the threshold is two molecules, or two
distant parts of one, occupying the same space.

**The threshold, and the evidence for it.**  On a healthy embedded example-17 membrane (136,514
atoms, 48,933 heavy) the closest non-bonded heavy-atom pair in the whole system was **further
than 2.2 A** -- the scan found nothing at all inside that radius.  Real chemistry does come
closer than 2.2 in other systems (a short strong hydrogen bond reaches ~2.4, metal coordination
~2.0), so the default sits at 1.8 A: below anything with a chemical reason to be there, far
above the sub-angstrom overlaps that actually diverge.  Raising it risks refusing a legitimate
structure; lowering it toward zero gives back the failure this exists to prevent.
"""
import logging

import numpy as np
from scipy.spatial import cKDTree

from .density_convergence import _read_coor_xyz

logger = logging.getLogger(__name__)

#: Heavy-atom separation (A) below which a non-bonded pair is an overlap rather than a contact.
DEFAULT_MIN_CONTACT = 1.8

#: Search radius (A) for candidate pairs.  Only pairs inside this are examined, so it bounds the
#: cost; it must exceed any threshold a caller may set.
_SEARCH_RADIUS = 3.0


def _bonded_exclusions(struct, n_atoms):
    """Adjacency list over atom INDICES, from the PSF bond list."""
    serial_to_i = {a.serial: i for i, a in enumerate(struct.atoms)}
    adj = [[] for _ in range(n_atoms)]
    for b in struct.bonds:
        i, j = serial_to_i.get(b.serial1), serial_to_i.get(b.serial2)
        if i is None or j is None:
            continue
        adj[i].append(j)
        adj[j].append(i)
    return adj


def worst_nonbonded_contact(struct, xyz, box=None, radius=_SEARCH_RADIUS):
    """Closest non-bonded heavy-atom pair.

    Parameters
    ----------
    struct : PSFContents
        Parsed with ``parse_topology=['bonds']``; the bonds are what make the 1-2/1-3 exclusion
        possible, and without them every bonded pair reads as an overlap.
    xyz : (natom, 3) array
    box : (3,) array or None
        Orthorhombic side lengths for the minimum-image convention; ``None`` for a
        non-periodic system.

    Returns
    -------
    (distance, i, j)
        ``(inf, None, None)`` when no non-bonded heavy pair lies within ``radius``.
    """
    atoms = struct.atoms
    if len(xyz) != len(atoms):
        raise ValueError(f'{len(xyz)} coordinates against {len(atoms)} atoms in the PSF')
    names = np.array([a.atomname for a in atoms], dtype=object)
    heavy = np.array([not str(n).startswith('H') for n in names])
    if not heavy.any():
        return float('inf'), None, None
    idx = np.where(heavy)[0]
    if box is not None:
        L = np.asarray(box, dtype=float)
        pts = np.mod(xyz[idx], L)
        tree = cKDTree(pts, boxsize=L)
    else:
        L, pts = None, xyz[idx]
        tree = cKDTree(pts)
    pairs = tree.query_pairs(r=radius, output_type='ndarray')
    if not len(pairs):
        return float('inf'), None, None
    adj = _bonded_exclusions(struct, len(atoms))
    worst = (float('inf'), None, None)
    for a, b in pairs:
        ia, ib = int(idx[a]), int(idx[b])
        if ib in adj[ia] or any(ib in adj[k] for k in adj[ia]):   # 1-2 or 1-3
            continue
        delta = pts[a] - pts[b]
        if L is not None:
            delta = (delta + L / 2) % L - L / 2
        d = float(np.linalg.norm(delta))
        if d < worst[0]:
            worst = (d, ia, ib)
    return worst


def describe_atom(struct, i):
    a = struct.atoms[i]
    return f'{a.segname} {a.resname} {a.resid} {a.atomname}'


def check_contacts(psf_path, coor_path, box=None, min_distance=DEFAULT_MIN_CONTACT):
    """``(ok, message)`` for the coordinates NAMD is about to read.

    ``coor_path`` must be the file NAMD will actually use -- the binary ``.coor`` when one
    exists, not the ``.pdb`` beside it, which may be an earlier stage's.  Checking the wrong
    file would be a gate that reports on coordinates nobody runs.
    """
    from ..psfutil.psfcontents import PSFContents
    struct = PSFContents(str(psf_path), parse_topology=['bonds'])
    xyz = _read_coor_xyz(coor_path, len(struct.atoms))
    d, i, j = worst_nonbonded_contact(struct, xyz, box=box)
    if d >= min_distance:
        shown = 'none within the search radius' if d == float('inf') else f'{d:.3f} A'
        return True, f'closest non-bonded heavy-atom pair: {shown}'
    return False, (f'{describe_atom(struct, i)} and {describe_atom(struct, j)} are {d:.3f} A '
                   f'apart (minimum {min_distance} A); they are not bonded to each other, so '
                   f'this is a steric overlap, and dynamics on it diverges -- usually as '
                   f'"Constraint failure in RATTLE algorithm". Minimize first, or rebuild')
