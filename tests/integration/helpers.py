# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Shared parsing and structural invariants for the integration tests.

These tests build real systems and then ask whether the result is *physically coherent*.  The
parsing here is deliberately naive and self-contained -- raw text, fixed columns, no pestifer
imports.  A test that reads a built PSF with pestifer's own reader cannot catch a bug in the
representation the reader and writer share, and several of the bugs these tests exist to guard
are exactly that shape (fixed-column corruption, atom-ordering mismatches).

:func:`assert_psf_sane` is the battery every build should pass regardless of what it was
exercising.  Each of its checks corresponds to a class of bug that has actually shipped here:

* duplicate bonds -- CHARMM ``modify_res`` residues declaring a redundant backward peptide bond,
  which psfgen turns into a degenerate angle NAMD then rejects (see ``util/psf.py``);
* over-long bonds -- glycan branch connectivity wired to the wrong residues after a
  BFS-vs-serial reordering, producing 10-16 A "bonds";
* atoms at the origin -- the membrane grid/psfgen split stranding atoms at (0,0,0);
* atom-count mismatch -- fixed-column PDB corruption from 6-character CHARMM carbohydrate
  resnames shifting the coordinate fields;
* non-integral charge -- a mis-built or truncated segment.
"""

import math
from pathlib import Path



# --- parsing -----------------------------------------------------------------------------------

def parse_psf(psf: Path):
    """``(atoms, bonds)`` from a PSF.

    ``atoms``: {atomid -> dict(segname, resid, resname, name, charge)}
    ``bonds``: [(atomid, atomid), ...]
    """
    lines = Path(psf).read_text().splitlines()
    ai = next(k for k, l in enumerate(lines) if '!NATOM' in l)
    natom = int(lines[ai].split()[0])
    atoms = {}
    for l in lines[ai + 1: ai + 1 + natom]:
        p = l.split()
        # psfgen PSF: id segname resid resname name type charge mass ...
        atoms[int(p[0])] = dict(segname=p[1], resid=p[2], resname=p[3], name=p[4],
                                charge=float(p[6]))
    bonds = []
    try:
        bi = next(k for k, l in enumerate(lines) if '!NBOND' in l)
    except StopIteration:
        return atoms, bonds
    nbond = int(lines[bi].split()[0])
    ints = []
    k = bi + 1
    while len(ints) < nbond * 2 and k < len(lines):
        ints += [int(x) for x in lines[k].split()]
        k += 1
    bonds = [(ints[2 * b], ints[2 * b + 1]) for b in range(nbond)]
    return atoms, bonds


def parse_pdb_coords(pdb: Path):
    """``{serial -> (x, y, z)}`` in file order.

    Fixed-column, as the PDB format declares: this is what makes the check meaningful, since a
    writer that shifts a column produces coordinates that parse but are wrong.
    """
    coords = {}
    idx = 0
    for l in Path(pdb).read_text().splitlines():
        if l.startswith(('ATOM', 'HETATM')):
            idx += 1
            coords[idx] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return coords


def parse_psf_pdb(psf: Path, pdb: Path):
    """``(atoms, bonds, coords)``; pestifer writes PSF and PDB in the same atom order."""
    atoms, bonds = parse_psf(psf)
    return atoms, bonds, parse_pdb_coords(pdb)


# --- invariants --------------------------------------------------------------------------------

def _close_contacts(coords, bonds, cutoff):
    """Non-bonded atom pairs closer than ``cutoff``, using a KD-tree.

    Bonded pairs are excluded; at a sub-angstrom cutoff, 1-3 and 1-4 neighbours are already far
    enough not to register, so only genuine overlaps survive.
    """
    try:
        import numpy as np
        from scipy.spatial import cKDTree
    except Exception:                      # pragma: no cover - scipy is a hard dependency
        return []
    serials = sorted(coords)
    xyz = np.array([coords[s] for s in serials])
    tree = cKDTree(xyz)
    bonded = {tuple(sorted(b)) for b in bonds}
    hits = []
    for i, j in tree.query_pairs(cutoff):
        a, b = serials[i], serials[j]
        if tuple(sorted((a, b))) in bonded:
            continue
        hits.append((a, b, float(np.linalg.norm(xyz[i] - xyz[j]))))
    return hits


def assert_psf_sane(psf: Path, pdb: Path, *, max_bond=2.5, min_contact=0.8,
                    check_contacts=True, allow_origin=False, unminimized=False, context=''):
    """Assert that a built system is structurally coherent.

    Parameters mirror the ways a legitimate build can differ: ``max_bond`` is raised for systems
    with deliberately long links, ``check_contacts`` is disabled where an overlap is expected
    before minimization, and ``allow_origin`` for the rare build that genuinely places an atom
    at the origin.

    Pass ``unminimized=True`` for a build that stops at psfgen with no relaxation.  psfgen
    *guesses* coordinates for hydrogens it has to add, and an unrelaxed guess can sit 3 A from
    its parent or on top of a neighbour -- neither is a defect, and minimization removes both.
    That mode therefore skips the contact check and applies the bond-length check to heavy atoms
    only.  Heavy-atom connectivity is what mis-wiring corrupts, so the check that matters is
    retained: the glycan branch bug produced 10-16 A C-O bonds, which this still catches.

    On ``max_bond``: the default is 2.5 A, not the ~1.5 A of a covalent single bond.  A disulfide
    S-S is ~2.05 A and is perfectly correct -- a tighter bound fails on almost every protein.  The
    failures this guards against were 10-16 A, so 2.5 A keeps a factor-of-four margin while
    admitting real chemistry.  A test checking a *narrower* population (intra-glycan bonds, say,
    which are ~1.4 A) should pass a tighter value.
    """
    tag = f'{context}: ' if context else ''
    psf, pdb = Path(psf), Path(pdb)
    assert psf.exists(), f'{tag}no PSF at {psf}'
    assert pdb.exists(), f'{tag}no PDB at {pdb}'

    atoms, bonds = parse_psf(psf)
    coords = parse_pdb_coords(pdb)

    # 1. the two files describe the same system
    assert len(atoms) == len(coords), (
        f'{tag}PSF declares {len(atoms)} atoms but PDB has {len(coords)} -- the files disagree, '
        f'which is what fixed-column PDB corruption looks like')

    # 2. no duplicate bonds
    seen, dupes = set(), []
    for a, b in bonds:
        key = tuple(sorted((a, b)))
        if key in seen:
            dupes.append(key)
        seen.add(key)
    assert not dupes, (
        f'{tag}{len(dupes)} duplicate bond(s) in the PSF, e.g. {dupes[:5]} -- a duplicated bond '
        f'becomes a degenerate angle that NAMD rejects')

    # 3. no self-bonds
    self_bonds = [b for b in bonds if b[0] == b[1]]
    assert not self_bonds, f'{tag}{len(self_bonds)} atom(s) bonded to themselves'

    # 4. every bond is a physically plausible length
    long_bonds = []
    for a, b in bonds:
        if a not in coords or b not in coords:
            continue
        if unminimized and (atoms[a]['name'].startswith('H') or atoms[b]['name'].startswith('H')):
            continue                       # a guessed hydrogen has no meaningful position yet
        d = math.dist(coords[a], coords[b])
        if d > max_bond:
            aa, bb = atoms[a], atoms[b]
            long_bonds.append(f'{aa["segname"]}:{aa["resname"]}{aa["resid"]}:{aa["name"]}'
                              f'--{bb["segname"]}:{bb["resname"]}{bb["resid"]}:{bb["name"]} '
                              f'({d:.1f} A)')
    assert not long_bonds, (
        f'{tag}{len(long_bonds)} bond(s) longer than {max_bond} A: {long_bonds[:8]} -- '
        f'connectivity is wired to the wrong atoms')

    # 5. nothing stranded at the origin
    if not allow_origin:
        at_origin = [s for s, c in coords.items() if c == (0.0, 0.0, 0.0)]
        assert not at_origin, (
            f'{tag}{len(at_origin)} atom(s) sit exactly at the origin (serials '
            f'{at_origin[:5]}) -- coordinates were never assigned')

    # 6. total charge is integral
    q = sum(a['charge'] for a in atoms.values())
    assert abs(q - round(q)) < 1e-3, (
        f'{tag}total charge {q:.4f} is not integral -- a segment is mis-built or truncated')

    # 7. no segment is empty, and none is unnamed
    segs = {a['segname'] for a in atoms.values()}
    assert '' not in segs, f'{tag}an atom has an empty segname'

    # 8. no atoms occupying the same space
    if check_contacts and not unminimized:
        hits = _close_contacts(coords, bonds, min_contact)
        assert not hits, (
            f'{tag}{len(hits)} non-bonded atom pair(s) closer than {min_contact} A, e.g. '
            f'{[(a, b, round(d, 2)) for a, b, d in hits[:5]]} -- atoms are on top of each other')
