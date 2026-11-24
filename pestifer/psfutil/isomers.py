# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Double-bond isomers: hold every C=C at its CHARMM isomer through an md stage, and stop a build
that ends a stage with one wrong.

Found 2026-11-24 on example 17, after the athermal-MC sampler and the ring resolver had both been
stopped from turning about a C=C (:mod:`pestifer.util.bondrules`).  The rebuilt quilts still
carried 13, 7 and 6 trans oleoyl bonds.  Traced through one replica's lower-leaflet patch, every
C=C was cis through psfgen, solvation and the ring check, and two were trans after the FIRST
MINIMIZE: the raw grid placement put other lipids' carbons 0.76-1.4 A from them (123 inter-lipid
heavy-atom pairs under 1 A in that patch), and the minimizer relieved the overlap by turning the
double bond.  ``X CEL1 CEL1 X`` has multiplicity 2, so trans is a minimum too, and no later MD
restores it.  ``write_grid_pdb`` warned about the overlaps and said "relaxation should resolve
these"; it did, by flipping the bond.

So this is the same failure the protein stereo hold (:mod:`pestifer.psfutil.stereo`) answers --
a minimizer started from a crowded placement does not care which way it relieves the strain --
and it gets the same two answers: hold the configuration through the stages that start from a
crowded structure, and check it afterwards.  The difference is the check: a wrong protein center
is reported, a wrong double bond stops the build, because nothing downstream (relaxation,
pressure profiles, a published APL) can tell a trans oleoyl from a cis one.

**The reference is CHARMM's, not an earlier stage.**  Each double bond's isomer is read from the
residue's internal coordinates -- a proper IC whose middle pair is the double bond, e.g. POPC's
``IC C28 C29 C210 C211 ... -1.49`` (cis) or PSM's ``IC C6S C5S C4S C3S ... 179.01`` (trans) --
which is what psfgen builds from.  A bond whose residue is not provisioned, or has no such IC, is
held at whatever isomer it has when the stage starts and is not checked.

Which bonds are double bonds is decided by the shared topological rule
(:func:`pestifer.util.bondrules.is_multiple_bond`, full coordination from the PSF).  Ring bonds
are left out: an aromatic ring satisfies the rule at every bond, is held planar by its own terms,
and has no isomer to lose.
"""
from __future__ import annotations

import logging
from collections import defaultdict, deque
from dataclasses import dataclass

import numpy as np

from ..util.bondrules import element_of_mass, is_multiple_bond

logger = logging.getLogger(__name__)

RING_SEARCH_DEPTH = 8
"""A C=C whose ends are rejoined by a path of at most this many bonds is a ring bond and skipped."""


@dataclass(frozen=True)
class IsomerDefect:
    segname: str
    resid: str
    resname: str
    names: tuple            # the four atom names of the dihedral
    reference: str          # 'cis' or 'trans', from the CHARMM IC
    value: float            # the dihedral now, degrees

    def describe(self):
        found = 'trans' if self.reference == 'cis' else 'cis'
        return (f'{self.segname}:{self.resname}{self.resid} {self.names[1]}={self.names[2]} is {found} '
                f'({"-".join(self.names)} {self.value:.0f} deg); CHARMM defines it {self.reference}')


def isomer_of(deg: float) -> str:
    return 'cis' if abs(deg) < 90.0 else 'trans'


def isomer_reference(charmmff_content, resnames) -> dict:
    """``{resname: {frozenset(middle pair): (n1, n2, n3, n4, degrees)}}`` from the residues' proper
    ICs.  A provisioned residue is taken as provisioned; any other is read from its topology file
    with :meth:`~pestifer.charmmff.charmmffcontent.CHARMMFFContent.read_resi`, never by
    provisioning, which would narrow the content another task relies on.  Found necessary the first
    time this ran end to end: a continuation build provisions nothing, and the check then had a
    reference for none of 114 double bonds and reported them all as right.  A name no topology file
    defines simply has no reference."""
    ref = {}
    if charmmff_content is None:
        return ref
    for rn in sorted(set(resnames)):
        resi = None
        try:
            if getattr(charmmff_content, 'provisioned', False):
                resi = charmmff_content.get_resi(rn)
            if resi is None:
                resi = charmmff_content.read_resi(rn)
        except Exception as e:
            logger.debug(f'{rn}: no CHARMM definition read for the isomer reference ({type(e).__name__}: {e})')
        if resi is None:
            continue
        table = {}
        for ic in getattr(resi, 'ICs', []):
            if ic.dihedral_type != 'proper':
                continue
            names = tuple(ic.atoms)
            if any(n[:1] in '+-*' or n[:1].isdigit() for n in names):
                continue
            table[frozenset(names[1:3])] = (*names, float(ic.dihedral.degrees))
        if table:
            ref[rn] = table
    return ref


def _read_psf(psf):
    """Lean PSF read: per-atom (segname, resid, resname, name, mass) and the bond list (0-based)."""
    with open(psf) as f:
        lines = f.read().split('\n')
    i = next(k for k, l in enumerate(lines) if '!NATOM' in l)
    n = int(lines[i].split()[0])
    atoms = []
    for l in lines[i + 1:i + 1 + n]:
        t = l.split()
        atoms.append((t[1], t[2], t[3], t[4], float(t[7])))
    j = next(k for k, l in enumerate(lines) if '!NBOND' in l)
    nb = int(lines[j].split()[0])
    vals = []
    k = j + 1
    while len(vals) < 2 * nb:
        vals.extend(int(x) - 1 for x in lines[k].split())
        k += 1
    return atoms, np.array(vals, dtype=int).reshape(-1, 2) if nb else np.zeros((0, 2), dtype=int)


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - np.sum(b0 * b1, axis=-1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, axis=-1, keepdims=True) * b1
    return np.degrees(np.arctan2(np.sum(np.cross(b1, v) * w, axis=-1), np.sum(v * w, axis=-1)))


class DoubleBondChecker:
    """Every non-ring C=C of one topology, each with the four atoms that measure its isomer and,
    where CHARMM gives one, the isomer it should have.  Parse once; :meth:`check_coords` and
    :meth:`restraint_terms` then cost an array read per frame.

    ``reference`` is :func:`isomer_reference`'s table, or a callable taking the set of resnames
    that carry a C=C and returning one, or None for no checking at all."""

    def __init__(self, psf, reference: dict | None = None):
        atoms, bonds = _read_psf(psf)
        self.natoms = len(atoms)
        nbr = defaultdict(list)
        for a, b in bonds:
            nbr[a].append(b)
            nbr[b].append(a)
        elem = [element_of_mass(a[4]) for a in atoms]
        self.terms = []       # (rows, segname, resid, resname, names, reference isomer or None)
        if not any(e == 'H' for e in elem):
            return            # no explicit hydrogens: coordination says nothing about bond order
        candidates = [(a, b) for a, b in bonds
                      if is_multiple_bond(elem[a], len(nbr[a]), elem[b], len(nbr[b]))
                      and not self._in_ring(nbr, a, b)]
        if callable(reference):
            reference = reference({atoms[a][2] for a, _ in candidates})
        reference = reference or {}
        rows_of = {}          # (seg, resid) -> {name: row}, built lazily per residue
        for a, b in candidates:
            seg, resid, resname = atoms[a][:3]
            key = (seg, resid)
            if key not in rows_of:
                rows_of[key] = {}
            rows = rows_of[key]
            if not rows:
                # every atom of this residue: walk the bonded molecule from a, keep same-residue atoms
                seen, todo = {a}, deque([a])
                while todo:
                    x = todo.popleft()
                    if atoms[x][:2] == key:
                        rows[atoms[x][3]] = x
                        for y in nbr[x]:
                            if y not in seen:
                                seen.add(y)
                                todo.append(y)
            quartet, ref = None, None
            ic = reference.get(resname, {}).get(frozenset((atoms[a][3], atoms[b][3])))
            if ic is not None and all(n in rows for n in ic[:4]) \
                    and atoms[b][:2] == key:
                quartet = tuple(rows[n] for n in ic[:4])
                ref = isomer_of(ic[4])
            if quartet is None:
                end_a = self._substituent(nbr, elem, a, b)
                end_b = self._substituent(nbr, elem, b, a)
                if end_a is None or end_b is None:
                    continue
                quartet = (end_a, a, b, end_b)
            self.terms.append((quartet, seg, resid, resname,
                               tuple(atoms[r][3] for r in quartet), ref))

    @staticmethod
    def _in_ring(nbr, a, b):
        depth = {a: 0}
        todo = deque([a])
        while todo:
            x = todo.popleft()
            if depth[x] >= RING_SEARCH_DEPTH:
                continue
            for y in nbr[x]:
                if x == a and y == b:
                    continue
                if y == b:
                    return True
                if y not in depth:
                    depth[y] = depth[x] + 1
                    todo.append(y)
        return False

    @staticmethod
    def _substituent(nbr, elem, x, partner):
        others = sorted(y for y in nbr[x] if y != partner)
        for pick in (lambda y: elem[y] == 'C', lambda y: elem[y] != 'H', lambda y: True):
            hit = [y for y in others if pick(y)]
            if hit:
                return hit[0]
        return None

    @property
    def nreferenced(self):
        return sum(1 for t in self.terms if t[5] is not None)

    def load_coords(self, path):
        """Coordinates from a PDB or a NAMD binary ``.coor``."""
        path = str(path)
        if path.endswith('.coor'):
            raw = open(path, 'rb').read()
            n = int(np.frombuffer(raw[:4], dtype='<i4')[0])
            xyz = np.frombuffer(raw[4:4 + 24 * n], dtype='<f8').reshape(n, 3).copy()
        else:
            xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                            for l in open(path) if l.startswith(('ATOM', 'HETATM'))])
        if xyz.shape[0] != self.natoms:
            raise ValueError(f'{path} has {xyz.shape[0]} atoms; the PSF has {self.natoms}')
        return xyz

    def angles(self, xyz):
        if not self.terms:
            return np.zeros(0)
        r = np.array([t[0] for t in self.terms])
        return _dihedral(xyz[r[:, 0]], xyz[r[:, 1]], xyz[r[:, 2]], xyz[r[:, 3]])

    def check_coords(self, xyz) -> list[IsomerDefect]:
        """Every C=C whose isomer differs from CHARMM's.  Bonds with no reference are not checked."""
        out = []
        for t, d in zip(self.terms, self.angles(xyz)):
            if t[5] is not None and isomer_of(d) != t[5]:
                out.append(IsomerDefect(t[1], t[2], t[3], t[4], t[5], float(d)))
        return out

    def restraint_terms(self, xyz):
        """``(rows, held, segname, resid, resname, names)`` per held bond: a referenced bond that is
        right now is held at its reference (0 or 180 deg); an unreferenced one at the isomer it has.
        A referenced bond that is already wrong is left out -- holding it would cement it, and
        pulling it back through the barrier is a repair, not a hold; the check reports it."""
        held = []
        for t, d in zip(self.terms, self.angles(xyz)):
            now = isomer_of(d)
            if t[5] is not None and now != t[5]:
                continue
            held.append((t[0], 0.0 if now == 'cis' else 180.0, t[1], t[2], t[3], t[4]))
        return held

    @staticmethod
    def term_strain(xyz, terms):
        """Signed deviation (degrees) of each held term from its held angle at ``xyz``."""
        if not terms:
            return np.zeros(0)
        r = np.array([t[0] for t in terms])
        ang = _dihedral(xyz[r[:, 0]], xyz[r[:, 1]], xyz[r[:, 2]], xyz[r[:, 3]])
        return (ang - np.array([t[1] for t in terms]) + 180.0) % 360.0 - 180.0
