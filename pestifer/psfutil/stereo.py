# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Protein stereochemistry check: D alpha carbons, inverted Thr/Ile beta carbons, and cis peptide
bonds that do not precede a proline.

None of these is something a build step ever asks for, and none of them is visible to anything
else pestifer checks.  Found 2026-10-29 from an outside group's audit of HIV Env systems: the
residues pestifer had modeled onto 4ZMJ's gp41 C terminus came out D-Trp666 and D-Lys683 in all
three protomers.  Example 11 shows the same thing on a current build: psfgen's output is clean,
and the first minimization makes the modeled B:Leu149 D, with a cis Asp148-Leu149 bond that
survives to the final structure.  A minimizer started from a crowded placement pushes an alpha
hydrogen through the N-C-CB plane as readily as it relieves the clash, and nothing pushes it
back.

Chirality is read as the signed volume of three substituents about the stereocenter:

    V = (N - CA) . ((C - CA) x (CB - CA))      L-amino acid: V > 0

and for the beta carbons of threonine (OG1, CG2, CA about CB) and isoleucine (CG1, CG2, CA about
CB), whose natural configurations also give V > 0.  The signs were read off a real build (all
1549 alpha carbons of example 11's psfgen output, and every Thr and Ile in it) rather than
derived, because a sign convention is exactly the kind of claim that is easy to get backwards.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np

from ..util.coord import coorddf_from_pdb
from .psfcontents import PSFContents

logger = logging.getLogger(__name__)

CIS_CUTOFF = 30.0
"""A peptide bond with ``|omega|`` below this many degrees is reported as cis."""

# stereocenter, then the three substituents in the order that gives a positive volume
_CENTERS = {
    'D-CA': ('CA', ('N', 'C', 'CB')),
    'THR-CB': ('CB', ('OG1', 'CG2', 'CA')),
    'ILE-CB': ('CB', ('CG1', 'CG2', 'CA')),
}
_APPLIES = {'D-CA': None, 'THR-CB': {'THR'}, 'ILE-CB': {'ILE'}}


@dataclass(frozen=True)
class StereoDefect:
    kind: str                  # 'D-CA', 'THR-CB', 'ILE-CB' or 'cis'
    segname: str
    resid: str
    resname: str
    value: float               # signed volume (A^3) or omega (degrees)

    @property
    def key(self):
        """Identity across frames of one topology: the value moves, the defect does not."""
        return (self.kind, self.segname, self.resid)

    def describe(self):
        what = {'D-CA': 'D alpha carbon', 'THR-CB': 'inverted Thr beta carbon (allo-Thr)',
                'ILE-CB': 'inverted Ile beta carbon (allo-Ile)'}.get(self.kind)
        if what is None:
            what = f'cis peptide bond to the preceding residue (omega {self.value:.0f} deg)'
        return f'{self.segname}:{self.resname}{self.resid} {what}'


def _dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - np.sum(b0 * b1, axis=-1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, axis=-1, keepdims=True) * b1
    x = np.sum(v * w, axis=-1)
    y = np.sum(np.cross(b1, v) * w, axis=-1)
    return np.degrees(np.arctan2(y, x))


class StereoChecker:
    """Prepared check for one topology; :meth:`check` reads one coordinate frame.

    The PSF is parsed once, so repeated frames of the same system cost only an array read.
    Residues missing any atom a test needs are skipped by that test, not guessed at.
    """

    def __init__(self, psf):
        self.topol = PSFContents(psf, parse_topology=['bonds'], topology_segtypes=['protein'])
        atoms = list(self.topol.atoms.data)
        self.natoms = len(atoms)
        row_of = {a.serial: i for i, a in enumerate(atoms)}
        residues = {}
        for i, a in enumerate(atoms):
            if a.segtype != 'protein':
                continue
            key = (a.segname, str(a.resid.resid))
            residues.setdefault(key, {'resname': a.resname, 'rows': {}})['rows'][a.atomname] = i
        self.centers = []
        for kind, (center, subs) in _CENTERS.items():
            only = _APPLIES[kind]
            for (seg, resid), r in residues.items():
                if only is not None and r['resname'] not in only:
                    continue
                if kind == 'D-CA' and r['resname'] == 'GLY':
                    continue
                rows = r['rows']
                if center in rows and all(s in rows for s in subs):
                    self.centers.append((kind, seg, resid, r['resname'],
                                         rows[center], [rows[s] for s in subs]))
        # peptide bonds: a PSF bond C(i)-N(j) between two protein residues of one segment
        self.peptides = []
        for b in self.topol.bonds:
            i, j = (row_of[s] for s in b.idx_list)
            ai, aj = atoms[i], atoms[j]
            if ai.segtype != 'protein' or aj.segtype != 'protein':
                continue
            if (ai.atomname, aj.atomname) == ('N', 'C'):
                i, j, ai, aj = j, i, aj, ai
            if (ai.atomname, aj.atomname) != ('C', 'N'):
                continue
            ri = residues.get((ai.segname, str(ai.resid.resid)), {}).get('rows', {})
            rj = residues.get((aj.segname, str(aj.resid.resid)), {}).get('rows', {})
            if 'CA' in ri and 'CA' in rj:
                self.peptides.append((aj.segname, str(aj.resid.resid), aj.resname,
                                      ri['CA'], i, j, rj['CA']))

    def load_coords(self, pdb):
        df = coorddf_from_pdb(pdb)
        if df.shape[0] != self.natoms:
            raise ValueError(f'{pdb} has {df.shape[0]} atoms; the PSF has {self.natoms}')
        return np.array(df[['x', 'y', 'z']].values, dtype=float)

    def check_coords(self, xyz) -> list[StereoDefect]:
        out = []
        if self.centers:
            c = xyz[[t[4] for t in self.centers]]
            s = np.array([t[5] for t in self.centers])
            a, b, d = (xyz[s[:, k]] - c for k in range(3))
            vol = np.sum(a * np.cross(b, d), axis=1)
            for t, v in zip(self.centers, vol):
                if v < 0:
                    out.append(StereoDefect(t[0], t[1], t[2], t[3], float(v)))
        if self.peptides:
            p = np.array([t[3:] for t in self.peptides])
            omega = _dihedral(xyz[p[:, 0]], xyz[p[:, 1]], xyz[p[:, 2]], xyz[p[:, 3]])
            for t, w in zip(self.peptides, omega):
                if abs(w) < CIS_CUTOFF and t[2] != 'PRO':
                    out.append(StereoDefect('cis', t[0], t[1], t[2], float(w)))
        return out

    def restraint_lines(self, xyz, k=100.0):
        """NAMD ``extraBonds`` lines holding every L alpha carbon L and every trans peptide bond trans
        at coordinates ``xyz``, with harmonic constant ``k``: the N-C-CB-CA improper at -33 deg (an L
        center reads about -33) and omega at 180 deg.  Row indices are 0-based, as NAMD wants.

        Centers or bonds already wrong are left alone: a restraint pulling a D center back through
        the plane is a different operation from holding a correct one, and is the tracker's to report.
        Measured on the 4ZMJ MPER build that went D at 666/683 (2026-10-29): 3711 restraints at k 100
        kept every center L through the first minimize."""
        lines = []
        for t in self.centers:
            if t[0] != 'D-CA':
                continue
            ca = t[4]
            n, c, cb = t[5]
            if float(_dihedral(xyz[[n]], xyz[[c]], xyz[[cb]], xyz[[ca]])[0]) < 0:
                lines.append(f'dihedral {n} {c} {cb} {ca} {k} -33.0')
        for t in self.peptides:
            ca_i, c_i, n_j, ca_j = t[3:]
            w = float(_dihedral(xyz[[ca_i]], xyz[[c_i]], xyz[[n_j]], xyz[[ca_j]])[0])
            if abs(w) > 90:
                lines.append(f'dihedral {ca_i} {c_i} {n_j} {ca_j} {k} 180.0')
        return lines

    def check(self, pdb) -> list[StereoDefect]:
        return self.check_coords(self.load_coords(pdb))


def stereo_check(psf, pdb) -> list[StereoDefect]:
    """Convenience wrapper: every stereochemistry defect in one frame of one system."""
    return StereoChecker(psf).check(pdb)
