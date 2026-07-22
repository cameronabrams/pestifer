# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Apply post-psfgen rigid-body / alignment coordinate operations to a built psf+pdb in pure numpy --
the replacement for the VMD/Tcl emitters ``write_rottrans`` / ``write_align`` /
``write_transfer_coords`` in :mod:`pestifer.scripters.vmd` (coordinate-port plan, Group B / Phase 2).

A built system's coordinates live only in its psf/pdb fileset, so each op loads the pdb into an
:class:`~pestifer.molecule.atom.AtomList` (identity + coords), pulls per-atom mass and bond
connectivity from the PSF, applies the transform with numpy, and writes the pdb back -- the
"load / transform / save per task" model.  The principal-axis ``orient`` directive is intentionally
left on the VMD path (it is unused in shipping configs and its eigensolver order/sign convention
resists faithful reproduction).

Selection support is the small closed subset of VMD atomselect syntax actually used by these
directives: ``all``, ``protein``, ``backbone``, ``name X``, ``chain X``, ``segid X``, ``resid N``,
``resid N to M``, ``fragment N``, joined only by ``and``.
"""
import logging
import re

import numpy as np

from pidibble.pdbparse import PDBParser

from .atom import AtomList
from .transform import Transform
from ..core.labels import Labels
from ..psfutil.loop_ccd import dihedral_deg
from ..psfutil.psfcontents import PSFContents
from ..util.coord import rotate_points_about_axis

logger = logging.getLogger(__name__)

#: the protein-backbone atom names of VMD's ``backbone`` singleword (protein part) -- the four
#: main-chain atoms plus the C-terminal carboxyl oxygens (CHARMM ``OT1``/``OT2``, PDB ``OXT``),
#: which VMD flags as backbone
_PROTEIN_BACKBONE = frozenset({'N', 'CA', 'C', 'O', 'OT1', 'OT2', 'OXT'})

#: unit vectors for the x/y/z axis keywords
_AXIS_VEC = {'x': np.array([1.0, 0.0, 0.0]),
             'y': np.array([0.0, 1.0, 0.0]),
             'z': np.array([0.0, 0.0, 1.0])}


class CoordManipulateError(Exception):
    """Raised when a coordinate op's selections are incongruent, a selection term is unsupported,
    or a rigid-body fragment is not fully disconnected -- the numpy analogue of the Tcl
    ``PESTIFER-ERROR`` guards."""


class CoordManipulator:
    """
    A built psf+pdb loaded for numpy coordinate manipulation.

    Parameters
    ----------
    psf : str or None
        Path to the PSF (supplies per-atom mass, segid, segtype, and bonds).  ``None`` loads a
        coordinate-only reference from the PDB alone (masses default to 1; ``segid`` falls back to
        the PDB chain, and ``protein``/``backbone`` use resname classification; ``fragment`` and the
        disconnection guard are unavailable without bonds).
    pdb : str
        Path to the PDB (supplies coordinates and per-atom chain identity).
    """

    def __init__(self, psf, pdb):
        self.atoms = AtomList.from_pdb(PDBParser(filepath=pdb).parse().parsed)
        n = len(self.atoms)
        if psf:
            psfc = PSFContents(psf, parse_topology=['bonds'])
            if len(psfc.atoms) != n:
                raise CoordManipulateError(
                    f'psf/pdb atom-count mismatch: {len(psfc.atoms)} (psf) vs {n} (pdb)')
            self.atoms.apply_psf_attributes(psfc.atoms)     # segname/serial/resid/resname
            self._psf = psfc
            self._mass = np.array([a.atomicwt for a in psfc.atoms.data], dtype=float)
            self._segtype = [a.segtype for a in psfc.atoms.data]
            self._name = [a.atomname for a in psfc.atoms.data]
        else:
            self._psf = None
            self._mass = np.ones(n, dtype=float)
            self._segtype = [Labels.segtype_of_resname.get(a.resname, '') for a in self.atoms.data]
            self._name = [a.name for a in self.atoms.data]
        self._fragment = None

    # ---- coordinate access ------------------------------------------------

    @property
    def coords(self) -> np.ndarray:
        return self.atoms.coords

    @coords.setter
    def coords(self, arr: np.ndarray):
        self.atoms.coords = arr

    def write_pdb(self, path: str):
        """Write the (transformed) coordinates back to a PDB, standard columns (coords pinned at
        31-54), matching what VMD's ``writepdb`` produced and what NAMD/psfgen read next."""
        self.atoms.write_pdb(path, dialect='standard')

    # ---- selection --------------------------------------------------------

    def select(self, sel: str) -> np.ndarray:
        """Translate a (restricted) VMD atomselection string into a boolean mask over the atoms."""
        sel = (sel or 'all').strip()
        mask = np.ones(len(self.atoms), dtype=bool)
        for term in re.split(r'\s+and\s+', sel):
            mask &= self._term_mask(term.strip())
        return mask

    def _term_mask(self, term: str) -> np.ndarray:
        atoms = self.atoms.data
        toks = term.split()
        if term == 'all':
            return np.ones(len(atoms), dtype=bool)
        if term == 'protein':
            return np.array([st == 'protein' for st in self._segtype])
        if term == 'backbone':
            return np.array([st == 'protein' and nm in _PROTEIN_BACKBONE
                             for st, nm in zip(self._segtype, self._name)])
        if len(toks) == 2 and toks[0] == 'name':
            return np.array([a.name == toks[1] for a in atoms])
        if len(toks) == 2 and toks[0] == 'chain':
            return np.array([a.chainID == toks[1] for a in atoms])
        if len(toks) == 2 and toks[0] == 'segid':
            return np.array([a.segname == toks[1] for a in atoms])
        if len(toks) == 2 and toks[0] == 'resid':
            n = int(toks[1])
            return np.array([a.resid.resseqnum == n for a in atoms])
        if len(toks) == 4 and toks[0] == 'resid' and toks[2] == 'to':
            lo, hi = int(toks[1]), int(toks[3])
            return np.array([lo <= a.resid.resseqnum <= hi for a in atoms])
        if len(toks) == 2 and toks[0] == 'fragment':
            return self._fragments() == int(toks[1])
        raise CoordManipulateError(f'unsupported atom-selection term: {term!r}')

    def _fragments(self) -> np.ndarray:
        """Per-atom VMD fragment id: connected components over PSF bonds, numbered in ascending
        first-atom (index) order -- VMD's fragment-numbering convention."""
        if self._psf is None:
            raise CoordManipulateError('fragment selection requires a PSF (bond connectivity)')
        if self._fragment is None:
            n = len(self.atoms)
            parent = list(range(n))

            def find(i):
                while parent[i] != i:
                    parent[i] = parent[parent[i]]
                    i = parent[i]
                return i

            for b in self._psf.bonds.data:
                ri, rj = find(b.serial1 - 1), find(b.serial2 - 1)
                if ri != rj:
                    parent[max(ri, rj)] = min(ri, rj)
            frag = np.empty(n, dtype=int)
            roots, nextid = {}, 0
            for i in range(n):
                r = find(i)
                if r not in roots:
                    roots[r], nextid = nextid, nextid + 1
                frag[i] = roots[r]
            self._fragment = frag
        return self._fragment

    # ---- geometry helpers -------------------------------------------------

    def _center(self, mask: np.ndarray, weighted: bool) -> np.ndarray:
        """Center of the masked atoms; mass-weighted (``measure center ... weight mass``) or
        geometric (``measure center``).  VMD uses ``abs`` of the weights."""
        pts = self.coords[mask]
        if weighted:
            return np.average(pts, axis=0, weights=np.abs(self._mass[mask]))
        return pts.mean(axis=0)

    def _require_disconnected(self, mask: np.ndarray, sel: str):
        """Hard-error unless the masked fragment is fully disconnected from the rest of the system
        (no bond crosses the selection boundary), so a rigid-body move cannot deform the molecule."""
        if self._psf is None:
            raise CoordManipulateError(f'transrot selection "{sel}" needs a PSF to verify disconnection')
        in_sel = mask                     # index i (0-based) selected?  serial = i+1
        for b in self._psf.bonds.data:
            if in_sel[b.serial1 - 1] != in_sel[b.serial2 - 1]:
                raise CoordManipulateError(
                    f'transrot selection "{sel}" is not fully disconnected '
                    f'(a bond crosses the selection boundary); refusing to apply a rigid-body transform')

    # ---- operations -------------------------------------------------------

    def apply_rottrans(self, rt):
        """Apply a :class:`~pestifer.objs.rottrans.RotTrans` (TRANS/ROT/AXISANGLE/ALIGN)."""
        sel = rt.sel or 'all'
        mask = self.select(sel)
        if sel != 'all':
            self._require_disconnected(mask, sel)
        coords = self.coords
        if rt.movetype == 'TRANS':
            coords[mask] += np.array([rt.x, rt.y, rt.z], dtype=float)
        elif rt.movetype == 'ROT':
            com = self._center(mask, weighted=True)
            coords[mask] = rotate_points_about_axis(coords[mask], com, _AXIS_VEC[rt.axis], rt.angle)
        elif rt.movetype == 'AXISANGLE':
            selI, selJ, selK = rt.axis_atoms
            rI = self._center(self.select(selI), weighted=False)
            rJ = self._center(self.select(selJ), weighted=False)
            rK = self._center(self.select(selK), weighted=False)
            axis = np.cross(rI - rJ, rJ - rK)
            coords[mask] = rotate_points_about_axis(coords[mask], rJ, axis, rt.angle)
        elif rt.movetype == 'ALIGN':
            self._apply_align_vectors(rt, mask, coords)
        else:
            raise CoordManipulateError(f'unknown transrot movetype {rt.movetype!r}')
        self.coords = coords

    def _apply_align_vectors(self, rt, mask, coords):
        """The minimal (roll-free) rotation carrying ``source`` onto ``target``, about the masked
        fragment's mass-weighted center -- mirrors ``_write_align`` in vmd.py."""
        src = self._align_vector(rt.source)
        tgt = self._align_vector(rt.target)
        src = src / np.linalg.norm(src)
        tgt = tgt / np.linalg.norm(tgt)
        cross = np.cross(src, tgt)
        sin = np.linalg.norm(cross)
        cos = float(np.dot(src, tgt))
        com = self._center(mask, weighted=True)
        if sin < 1.0e-6:
            if cos < 0.0:                 # antiparallel: 180 about an arbitrary perpendicular
                perp = np.cross(src, [1.0, 0.0, 0.0])
                if np.linalg.norm(perp) < 1.0e-6:
                    perp = np.cross(src, [0.0, 1.0, 0.0])
                coords[mask] = rotate_points_about_axis(coords[mask], com, perp, 180.0)
            # else already aligned: no rotation
        else:
            angle = np.degrees(np.arctan2(sin, cos))
            coords[mask] = rotate_points_about_axis(coords[mask], com, cross, angle)

    def _align_vector(self, spec) -> np.ndarray:
        """An ALIGN vector: a literal ``[x, y, z]`` or a pair of selections ``[selA, selB]`` whose
        mass-weighted centers define the vector (from A to B)."""
        if len(spec) == 3 and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in spec):
            return np.array(spec, dtype=float)
        selA, selB = spec
        return self._center(self.select(selB), weighted=True) - self._center(self.select(selA), weighted=True)

    def apply_align(self, align, reference: "CoordManipulator"):
        """Least-squares-fit ``mobile_sel`` (this system) onto ``ref_sel`` (``reference``) and move
        ``apply_to`` by the fit -- the numpy ``measure fit`` / ``$mover move`` of :meth:`write_align`."""
        mob_mask = self.select(align.mobile_sel)
        ref_sel = align.ref_sel if align.ref_sel is not None else align.mobile_sel
        ref_mask = reference.select(ref_sel)
        nm, nr = int(mob_mask.sum()), int(ref_mask.sum())
        if nm != nr:
            raise CoordManipulateError(
                f'align: mobile_sel "{align.mobile_sel}" has {nm} atoms but ref_sel "{ref_sel}" has {nr}')
        transform, _ = Transform.superpose(self.coords[mob_mask], reference.coords[ref_mask])
        apply_mask = self.select(align.apply_to)
        coords = self.coords
        coords[apply_mask] = transform.apply(coords[apply_mask])
        self.coords = coords

    def apply_transfer_coords(self, tc, donor: "CoordManipulator"):
        """Copy ``donor_sel`` coordinates from ``donor`` onto ``mobile_sel`` here, optionally after
        rigidly pre-fitting the whole donor -- the numpy port of :meth:`write_transfer_coords`."""
        if tc.pre_align:
            d_mask = donor.select(tc.align_donor_sel)
            m_mask = self.select(tc.align_mobile_sel)
            nd, nm = int(d_mask.sum()), int(m_mask.sum())
            if nd != nm:
                raise CoordManipulateError(
                    f'transfer_coords align: align_donor_sel "{tc.align_donor_sel}" has {nd} atoms '
                    f'but align_mobile_sel "{tc.align_mobile_sel}" has {nm}')
            transform, _ = Transform.superpose(donor.coords[d_mask], self.coords[m_mask])
            donor.coords = transform.apply(donor.coords)     # move the entire donor
        d_sel = donor.select(tc.donor_sel)
        m_sel = self.select(tc.mobile_sel)
        nd, nm = int(d_sel.sum()), int(m_sel.sum())
        if nd != nm:
            raise CoordManipulateError(
                f'transfer_coords: donor_sel "{tc.donor_sel}" has {nd} atoms but '
                f'mobile_sel "{tc.mobile_sel}" has {nm}')
        coords = self.coords
        coords[m_sel] = donor.coords[d_sel]
        self.coords = coords

    def apply_orient(self, orient):
        """
        Reorient the whole molecule so its long principal axis lies along the target axis -- the
        numpy port of the VMD ``Orient`` package (mass-weighted inertia tensor -> principal axes).
        With a ``refatom`` the molecule is then recentered at the origin (geometric center) and
        flipped 180 about x when the reference atom sits at negative z.

        VMD's underlying eigensolver (``mevsvd_br``) returns the principal axes in an unsorted,
        sign-arbitrary order; this port instead aligns the least-inertia (long) axis with a
        deterministic sign, so the result is a canonical equivalent -- identical up to the roll
        about the target axis, which is physically immaterial for the membrane-embedding use this
        serves.
        """
        m = np.abs(self._mass)
        coords = self.coords
        com = np.average(coords, axis=0, weights=m)
        r = coords - com
        x, y, z = r[:, 0], r[:, 1], r[:, 2]
        inertia = np.array([
            [np.sum(m * (y * y + z * z)), -np.sum(m * x * y), -np.sum(m * x * z)],
            [-np.sum(m * x * y), np.sum(m * (x * x + z * z)), -np.sum(m * y * z)],
            [-np.sum(m * x * z), -np.sum(m * y * z), np.sum(m * (x * x + y * y))],
        ])
        _, evecs = np.linalg.eigh(inertia)          # ascending eigenvalues; columns are axes
        axisvec = evecs[:, 0]                        # smallest moment of inertia = molecular long axis
        target = _AXIS_VEC[orient.axis]
        if np.dot(axisvec, target) < 0.0:            # deterministic sign: rotate the short way
            axisvec = -axisvec
        cross = np.cross(axisvec, target)
        sin = np.linalg.norm(cross)
        if sin > 1.0e-9:
            angle = np.degrees(np.arctan2(sin, float(np.dot(axisvec, target))))
            coords = rotate_points_about_axis(coords, com, cross, angle)
        if orient.refatom:
            coords = coords - coords.mean(axis=0)    # recenter geometric center to the origin
            ref = self.select(f'name {orient.refatom}')
            if ref.any() and float(coords[ref][:, 2].mean()) < 0.0:
                coords = rotate_points_about_axis(coords, np.zeros(3), _AXIS_VEC['x'], 180.0)
        self.coords = coords

    # ---- internal-coordinate (torsion) rotations -- crot ------------------

    #: alpha-helix backbone target dihedrals (deg) used by fold_alpha
    _ALPHA = {'phi': -57.0, 'psi': -47.0, 'omega': 180.0}

    def _residue_index(self) -> np.ndarray:
        """Per-atom VMD ``residue`` index: a global 0-based counter incremented at each residue
        boundary (change of segname/resid/insertion) in atom order.  Verified to match VMD."""
        if getattr(self, '_ridx', None) is None:
            ridx = np.empty(len(self.atoms), dtype=int)
            r, prev = -1, None
            for i, a in enumerate(self.atoms.data):
                key = (a.segname, a.resid.resseqnum, a.resid.insertion)
                if key != prev:
                    r += 1
                    prev = key
                ridx[i] = r
            self._ridx = ridx
            self._names = np.array(self._name)
            # chain (pdb) of each residue index -- first atom of each residue
            nres = r + 1
            self._rchain = [None] * nres
            for i, a in enumerate(self.atoms.data):
                if self._rchain[ridx[i]] is None:
                    self._rchain[ridx[i]] = a.chainID
        return self._ridx

    def _residue_of(self, chain, resid) -> int:
        """VMD residue index of ``chain``/``resid`` (matched at its CA, else any atom)."""
        ridx = self._residue_index()
        cand = None
        for i, a in enumerate(self.atoms.data):
            if a.chainID == chain and a.resid.resseqnum == resid:
                if a.name == 'CA':
                    return int(ridx[i])
                if cand is None:
                    cand = int(ridx[i])
        if cand is None:
            raise CoordManipulateError(f'crot: residue chain {chain} resid {resid} not found')
        return cand

    def _same_chain(self, q, r) -> bool:
        ridx = self._residue_index()
        nres = len(self._rchain)
        if q < 0 or r < 0 or q >= nres or r >= nres:
            return False
        return self._rchain[q] == self._rchain[r]

    def _res_xyz(self, coords, ridx, names, r, atom_name):
        m = (ridx == r) & (names == atom_name)
        idx = np.nonzero(m)[0]
        return coords[idx[0]] if idx.size else None

    def apply_crot(self, crot, chainIDmap=None):
        """Apply a :class:`~pestifer.objs.crot.Crot` internal-coordinate rotation (PHI/PSI/OMEGA/
        CHI1/CHI2 via bond rotation, ALPHA via fold_alpha), remapping the chain through
        ``chainIDmap`` for biological-assembly images."""
        chain = (chainIDmap or {}).get(crot.chainID, crot.chainID)
        angle = crot.angle
        if angle in ('PHI', 'PSI', 'OMEGA'):
            r0 = self._residue_of(chain, crot.resid1.resseqnum)
            r1 = self._residue_of(chain, crot.resid2.resseqnum)
            direction = 'C' if crot.resid1 <= crot.resid2 else 'N'
            self._brot(r0, r1, angle.lower(), direction, crot.degrees)
        elif angle in ('CHI1', 'CHI2'):
            r0 = self._residue_of(chain, crot.resid1.resseqnum)
            self._brot(r0, -1, 'chi', int(angle[-1]), crot.degrees)
        elif angle == 'ALPHA':
            rb = self._residue_of(chain, crot.resid1.resseqnum)
            re = self._residue_of(chain, crot.resid2.resseqnum)
            rt = self._residue_of(chain, crot.resid3.resseqnum)
            self._fold_alpha(rb, re, rt)
        elif angle == 'ANGLEIJK':
            raise CoordManipulateError('ANGLEIJK is deprecated; use a transrot AXISANGLE instead')
        else:
            raise CoordManipulateError(f'crot angle {angle!r} is not supported on the numpy path')

    def _brot(self, r0, r1, angle_name, rot, deg):
        """Rotate the moving set about the backbone (or chi) bond by ``deg`` degrees -- the numpy
        port of crot.tcl ``brot``.  The moving set and rotation axis are chosen exactly as in the
        Tcl (residue-index ranges + main-chain-atom exclusions); rotation is about the first axis
        atom along the (second - first) direction (VMD ``trans bond``)."""
        ridx = self._residue_index()
        names = self._names
        coords = self.coords

        def xyz(r, nm):
            return self._res_xyz(coords, ridx, names, r, nm)

        if angle_name == 'phi':
            pn, pc = xyz(r0, 'N'), xyz(r0, 'CA')
            nter = np.isin(names, ['N', 'HN', 'HT1', 'HT2', 'HT3'])
            if rot == 'C':
                mask = ((ridx == r0) & ~nter) | ((ridx > r0) & (ridx <= r1))
            else:
                mask = ((ridx == r0) & nter) | ((ridx < r0) & (ridx >= r1))
        elif angle_name == 'psi':
            pn, pc = xyz(r0, 'CA'), xyz(r0, 'C')
            cter = np.isin(names, ['C', 'O', 'OT1', 'OT2', 'OXT'])
            if rot == 'C':
                mask = ((ridx == r0) & cter) | ((ridx > r0) & (ridx <= r1))
            else:
                mask = ((ridx == r0) & ~cter) | ((ridx < r0) & (ridx >= r1))
        elif angle_name == 'omega':
            pn, pc = xyz(r0, 'C'), xyz(r0 + 1, 'N')
            if rot == 'C':
                mask = (ridx > r0) & (ridx <= r1)
            else:
                mask = (ridx <= r0) & (ridx >= r1)
        elif angle_name == 'chi':
            if rot == 1:
                pn, pc = xyz(r0, 'CA'), xyz(r0, 'CB')
                mask = (ridx == r0) & ~np.isin(names, ['N', 'HN', 'CA', 'C', 'O'])
            else:
                pn, pc = xyz(r0, 'CB'), xyz(r0, 'CG')
                mask = (ridx == r0) & ~np.isin(names, ['N', 'HN', 'CA', 'CB', 'C', 'O'])
        else:
            raise CoordManipulateError(f'brot: angle {angle_name!r} not recognized')
        if pn is None or pc is None:
            raise CoordManipulateError(f'brot: missing axis atoms for {angle_name} at residue {r0}')
        coords[mask] = rotate_points_about_axis(coords[mask], pn, pc - pn, deg)
        self.coords = coords

    def residue_of_segname(self, segname, resid) -> int:
        """VMD residue index of ``segname``/``resid`` (matched at its CA, else any atom)."""
        ridx = self._residue_index()
        cand = None
        for i, a in enumerate(self.atoms.data):
            if a.segname == segname and a.resid.resseqnum == resid:
                if a.name == 'CA':
                    return int(ridx[i])
                if cand is None:
                    cand = int(ridx[i])
        if cand is None:
            raise CoordManipulateError(f'residue segname {segname} resid {resid} not found')
        return cand

    def apply_scrot(self, chi, residue_index, deg):
        """
        Rotate a side chain about its chi1 (CA-CB) or chi2 (CB-CG) bond -- the numpy port of
        crot.tcl ``SCrot_chi1``/``SCrot_chi2`` (rotation about the first axis atom, along the
        first-minus-second direction).  Unlike the ``brot`` CHI path this excludes HA/HB and CB
        from the moving set.
        """
        ridx = self._residue_index()
        names = self._names
        coords = self.coords

        def xyz(nm):
            idx = np.nonzero((ridx == residue_index) & (names == nm))[0]
            return coords[idx[0]] if idx.size else None

        if chi == 1:
            p1, p2 = xyz('CA'), xyz('CB')
            mask = (ridx == residue_index) & ~np.isin(names, ['N', 'HN', 'CA', 'CB', 'HA', 'C', 'O'])
        else:
            p1, p2 = xyz('CB'), xyz('CG')
            mask = (ridx == residue_index) & ~np.isin(names, ['N', 'HN', 'CA', 'CB', 'HA', 'C', 'O', 'HB1', 'HB2'])
        if p1 is None or p2 is None:
            raise CoordManipulateError(f'SCrot chi{chi}: missing axis atoms at residue {residue_index}')
        coords[mask] = rotate_points_about_axis(coords[mask], p1, p1 - p2, deg)
        self.coords = coords

    def _get_phi_psi_omega(self, r):
        ridx = self._residue_index()
        names = self._names
        coords = self.coords

        def xyz(rr, nm):
            return self._res_xyz(coords, ridx, names, rr, nm)

        phi = psi = omega = None
        if self._same_chain(r - 1, r):
            phi = dihedral_deg(xyz(r - 1, 'C'), xyz(r, 'N'), xyz(r, 'CA'), xyz(r, 'C'))
        if self._same_chain(r, r + 1):
            psi = dihedral_deg(xyz(r, 'N'), xyz(r, 'CA'), xyz(r, 'C'), xyz(r + 1, 'N'))
            omega = dihedral_deg(xyz(r, 'CA'), xyz(r, 'C'), xyz(r + 1, 'N'), xyz(r + 1, 'CA'))
        return phi, psi, omega

    def _fold_alpha(self, rbegin, rend, rterm):
        """Fold residues ``rbegin..rend`` toward alpha-helical phi/psi/omega by delta-rotating each
        current angle to its target -- the numpy port of crot.tcl ``fold_alpha``."""
        if rbegin < rend:
            side, inc = 'C', 1
            done = lambda r: r > rend
        else:
            side, inc = 'N', -1
            done = lambda r: r < rend
        r = rbegin
        while not done(r):
            phi, psi, omega = self._get_phi_psi_omega(r)
            for nm, cur in (('phi', phi), ('psi', psi), ('omega', omega)):
                if cur is None:
                    continue
                d = self._ALPHA[nm] - cur
                if side == 'N':
                    d = -d
                self._brot(r, rterm, nm, side, d)
            r += inc
