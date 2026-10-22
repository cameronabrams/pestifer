#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the Ring class and the ring_check function for detection of pierced rings
"""

import logging
import networkx as nx
import numpy as np
import time

from scipy.spatial import cKDTree

from functools import singledispatchmethod
from itertools import pairwise

from ..psfutil.psfbond import PSFBondList
from ..psfutil.psfcontents import PSFContents
from ..psfutil.psftopoelement import PSFTopoElementList,PSFTopoElement

from ..util.coord import coorddf_from_pdb, lawofcos
from ..util.linkcell import Linkcell
from ..util.util import countTime, cell_from_xsc

logger=logging.getLogger(__name__)

#: A single token carried by EVERY message that reports a ring piercing as a build failure.
#:
#: There are five distinct phrasings -- terminate's refusal, the ``delete: none`` stop, the
#: verify-pass survival, the protein/glycan-after-embed case -- and a consumer scraping build
#: logs had to know all of them.  On the 3.29.1 sweep it knew two: the classifier matched
#: nothing, a `grep` for findings reported **0 across 96 builds**, and the one real piercing
#: surfaced only through a looser fallback grep in the same readout.  A check that cannot
#: match is a check that silently did not run, which is this repository's recurring failure
#: wearing an operational costume.
#:
#: Grep for this, not for prose.  It is deliberately ugly so it cannot occur by accident and
#: cannot be reworded by someone improving a sentence.
#:
#: **And it contains no regex metacharacters, which is the second half of the lesson.**  The
#: first version was ``[RING-PIERCING]``.  That is a character CLASS matching any one of
#: ``R I N G - P E C``, so on a real build log a plain ``grep`` matched **88 of 94 lines** --
#: the banner, the author's email address, nearly everything.  Only ``grep -F`` was correct.
#:
#: That is the exact mirror of the failure it was introduced to fix: the old patterns could
#: only return zero, and the marker that replaced them could only return non-zero.  Neither
#: can fail, and a check that cannot fail is not a check.  A marker whose correctness depends
#: on the consumer remembering a flag is not a marker, so this one works under plain ``grep``,
#: ``grep -E``, ``grep -F``, ``awk`` and :mod:`re` alike.  ``test_the_marker_is_regex_safe``
#: pins that property rather than the spelling.
PIERCING_MARKER = 'RING_PIERCING_FOUND' 


class PSFRing(PSFTopoElement):
    """
    A class for handling rings in a molecular structure.
    This class represents a ring defined by a list of atom indices and provides methods to check if a bond pierces the ring.
    It also provides methods to yield treadmilled versions of the ring and to check for equality with another ring.

    Parameters
    ----------
    idx_list : list
        A list of atom indices that define the ring.
    """
    
    def treadmill(self):
        """ 
        yield the treadmilled versions of the list if atom indices in the ring
        This method generates all possible rotations of the ring's atom indices.
        Each rotation is a new arrangement of the atom indices, simulating the effect of a treadmill.
        
        Yields
        -------
        list
            A list of all possible rotations of the ring's atom indices.
        """
        for i in range(1,len(self.idx_list)):
            yield self.idx_list[i:]+self.idx_list[:i]

    def __eq__(self,other):
        """
        Check if this ring is equal to another ring.
        Two rings are considered equal if they have the same atom indices in any order, including their
        treadmilled versions.
        
        Parameters
        ----------
        other : Ring
            The other ring to compare with. 
        
        Returns
        -------
        bool
            True if the rings are equal, False otherwise.
        """
        check1=any([self.idx_list==other.idx_list] + [self.idx_list==x for x in other.treadmill()])
        check2=any([self.idx_list==other.idx_list[::-1]]+[self.idx_list==x[::-1] for x in other.treadmill()])
        return check1 or check2

    def pierced_by(self,bond,cutoff=3.5,tol=1.e-5):
        """
        Check if a bond pierces the ring.
        A bond is considered to pierce the ring if:

        - The center of mass (COM) of the bond is within a specified cutoff distance from the COM of the ring.
        - The two bond atoms lie on opposite sides of the ring's mean plane (the bond straddles it).
        - The ring's points, when projected onto a plane perpendicular to the bond at the bond's midpoint, form a closed loop (i.e., the sum of angles around the COM is approximately 2π).

        The second condition used to be a DISTANCE gate -- reject unless the bond midpoint lay
        within half a bond length of the ring COM -- under a reason string claiming it tested
        "both atoms on same side of ring plane".  It did not.  It rejected every threading that
        crossed the ring off-center, which is exactly what a small displacement of a pierced ring
        produces.  Measured 2026-10-22 on example 17's replica-1 quilt: the 3.29.4 resolver
        reported 6 piercings cleared, and 3 of them were still threaded -- each bond still
        straddling its ring plane with a winding of exactly 2π, but now crossing 1-2 A off-center,
        past a gate of 1.26 A.  Two had been "cleared" by sliding the sterol, one by spinning it.
        The resolver judges its own candidate poses with this function, so the gate let it
        certify its own non-fixes, and the terminate gate inherited the same blind spot.

        The straddle test is what the old reason string described.  Together with the winding
        test it means what it should: the bond SEGMENT crosses the ring plane INSIDE the ring.

        Parameters
        ----------
        bond : PSFBond
            The bond to check for piercing the ring.
        cutoff : float
            The cutoff distance for considering a bond to pierce the ring.
        tol : float
            The tolerance for checking if the ring forms a closed loop.

        Returns
        -------
        dict
            A dictionary containing the result of the piercing check. Keys include:

            - `pierced`: bool
                True if the bond pierces the ring, False otherwise.
            - `piercepoint`: np.ndarray
                The center of mass of the ring projected onto the plane perpendicular to the bond at the bond's midpoint, if the bond pierces the ring.
            - `error`: float
                The error in the closed loop condition, if the bond pierces the ring.
            - `reason`: str
                A reason for why the bond does not pierce the ring, if applicable.
        """
        dist=np.linalg.norm(bond.COM-self.COM)
        if dist>cutoff:
            return dict(pierced=False,reason='cutoff')
        # straddle: the two bond atoms must sit on opposite sides of the ring's mean plane.
        # Newell's normal (sum of cross products of consecutive ring radii) is robust for the
        # puckered rings this sees -- sugar chairs, sterol A/B rings -- and costs no SVD.
        R=self.P-self.COM
        normal=np.sum(np.cross(R,np.roll(R,-1,axis=0)),axis=0)
        nn=np.linalg.norm(normal)
        if nn==0.0:
            return dict(pierced=False,reason='degenerate ring')
        normal=normal/nn
        hA=float(np.dot(bond.P[0]-self.COM,normal))
        hB=float(np.dot(bond.P[1]-self.COM,normal))
        if hA*hB>=0.0:
            return dict(pierced=False,reason='both atoms on same side of ring plane')

        # project all ring points into plane perpendicular to bond at bond midpoint
        lbv,rbv=bond.radii[0],bond.radii[1]
        projected_ring=np.array([p+np.dot(lbv,bond.COM-p)*lbv for p in self.P])
        pCOM=np.mean(projected_ring,axis=0)
        phi=np.array([np.arccos(lawofcos(p[0]-bond.COM,p[1]-bond.COM)) for p in pairwise(np.concatenate((projected_ring,[projected_ring[0]])))])
        phisum=np.sum(phi)
        diff=phisum-np.pi*2
        if np.abs(diff)<tol:
            return dict(pierced=True,piercepoint=pCOM,error=np.abs(diff))
        return dict(pierced=False,reason='non-winding',error=np.abs(diff))

class RingList(PSFTopoElementList):
    """
    A class for handling lists of Ring objects.
    This class inherits from PSFTopoElementList and provides methods to initialize the list from a networkx graph,
    ingest coordinates, and check for string representation.
    
    Parameters
    ----------
    input_obj : nx.Graph or list of Ring
    """

    @singledispatchmethod
    def __init__(self,input_obj):
        self.data=input_obj
    
    @__init__.register(nx.Graph)
    def _from_graph(self,G,length_bound=None):
        # length_bound keeps the cycle search cheap and macrocycle-free: a protein
        # backbone closed by a disulfide is a huge chordless cycle, but every real
        # chemical ring (aromatic 5/6, proline 5, sugar 6, sterol 5/6) is <= 6, so
        # bounding to ~7 finds exactly those and skips the giant loops.
        L=[]
        for ll in nx.chordless_cycles(G,length_bound=length_bound):
            L.append(PSFRing(ll))
        super().__init__(L)

    def __str__(self):
        return ';'.join([str(x) for x in self])
    
class RingChecker:
    """Prepared pierced-ring checker.

    The coordinate-*independent* work -- parsing the PSF, building the bond list, and
    finding the (length-bounded) ring cycles -- is done once in the constructor, so a
    sequence of coordinate frames (e.g. trial side-chain rotamers in
    :class:`~pestifer.tasks.ringcheck.RingCheckTask`) can be checked with :meth:`check`
    without re-parsing the topology or re-finding the rings each time.

    Parameters
    ----------
    psf : str
        Path to the PSF file.
    cutoff : float
        Candidate search radius in Angstroms (bond midpoint to ring COM).  Any value at or
        above :meth:`PSFRing.pierced_by`'s 3.5 Å gate yields identical piercings; the default
        (4.0) sits just above it to prune candidates and keep the scan fast on big membranes.
    segtypes : list
        Segment types whose rings are checked.
    max_ring_size : int
        Only rings of at most this many atoms are considered.
    """

    def __init__(self, psf, cutoff=4.0, segtypes=['lipid'], max_ring_size=7):
        # cutoff is the candidate search radius; any value >= pierced_by's 3.5 A gate gives
        # identical results, so it stays just above the gate to keep the scan fast
        self.cutoff = cutoff
        self.topol = PSFContents(psf, parse_topology=['bonds'], topology_segtypes=segtypes)
        self.segname_to_segtype = {a.segname: a.segtype for a in self.topol.atoms}
        # key on str(resid): PSF atoms carry an int resid but the ring/bond resid comes
        # from the ingested coordinate frame as a string, so the tuples must be normalized
        self.resname_of = {(a.segname, str(a.resid)): a.resname for a in self.topol.atoms}
        self.rings = RingList(self.topol.G, length_bound=max_ring_size)
        # Precompute the coordinate-independent maps used by the fast targeted path
        # (:meth:`_check_fast`): atom-serial -> row, and per-ring / per-bond row indices,
        # so a targeted re-check does only vectorized array lookups (no per-element pandas
        # ``.loc`` and no link cell), turning a ~50 s whole-system scan into ~2 s.
        self._atoms = list(self.topol.atoms.data)
        self._serials = [a.serial for a in self._atoms]
        row_of = {s: i for i, s in enumerate(self._serials)}
        self._row_of_serial = row_of
        # heavy-atom mask (hydrogens named H...) for clash scoring
        self._heavy = np.array([not a.atomname.startswith('H') for a in self._atoms])
        self._bonds = list(self.topol.bonds)
        if self._bonds:
            self._bond_rows = np.array([[row_of[b.idx_list[0]], row_of[b.idx_list[1]]]
                                        for b in self._bonds], dtype=int)
            self._bond_meta = [(self._atoms[r[0]].segname, self._atoms[r[0]].resid.resid)
                               for r in self._bond_rows]
        else:
            self._bond_rows = np.empty((0, 2), dtype=int)
            self._bond_meta = []
        # cache each ring's row indices and derive its segname/resid from its first atom so
        # only_piercees filtering works without a full coordinate ingest
        self._ring_index = {}
        for ring in self.rings.data:
            ring._rows = np.array([row_of[s] for s in ring.idx_list], dtype=int)
            a0 = self._atoms[ring._rows[0]]
            ring.segname = a0.segname
            ring.resid = a0.resid.resid
            self._ring_index.setdefault((str(a0.segname), str(a0.resid.resid)), []).append(ring)
        logger.debug(f'RingChecker: parsed topology once, {len(self.rings)} rings')

    def check(self, pdb, xsc=None, only_piercees=None):
        """Check one coordinate frame for pierced rings and return the piercespecs.

        Only coordinates are re-read here; the PSF, bond list, and ring cycles from the
        constructor are reused.  ``xsc`` sets the periodic box; the system is treated as vacuum
        (box from coordinate extents) both when ``xsc`` is ``None`` and when the file carries no
        cell.  ``only_piercees`` (iterable of ``(segname, resid)``) restricts which rings are
        tested -- a cheap targeted re-check after a trial rotation.
        """
        topol = self.topol
        coorddf = coorddf_from_pdb(pdb)
        assert coorddf.shape[0] == len(topol.atoms), f'{pdb} is incongruent with the PSF'
        # guard on the parsed cell, not on the xsc path: cell_from_xsc returns (None, None)
        # for any xsc that carries no cell -- a vacuum run writes an origin-only
        # `step o_x o_y o_z`, 4 columns -- so a file that exists is not a promise of a box.
        box, orig = cell_from_xsc(xsc) if xsc is not None else (None, None)
        if only_piercees is not None:
            # targeted re-check after a trial rotation: vectorized, no link cell.  Only a
            # handful of named rings are tested, so skip the whole-system coordinate ingest.
            coords = coorddf[['x', 'y', 'z']].values
            return self._check_fast(coords, box, only_piercees)
        if box is not None:
            sidelengths = np.diagonal(box)
            ll = orig - 0.5 * sidelengths
            ur = orig + 0.5 * sidelengths
        else:
            ll = ur = None
            logger.debug('No periodic cell — treating system as non-periodic (vacuum)')
        coords = coorddf[['x', 'y', 'z']].values
        return self._scan(coords, box, self.rings.data, ll=ll, ur=ur)

    def _scan(self, coords, box, rings_to_check, ll=None, ur=None):
        """Whole-system pierced-ring scan over ``rings_to_check`` using a link cell for
        spatial acceleration.  Coordinates are taken from the ``coords`` array (in PSF atom
        order) -- bond and ring positions are set by vectorized numpy indexing rather than a
        per-element pandas ``.loc`` ingest, so a several-hundred-thousand-atom membrane scans
        in seconds instead of minutes.  ``ll``/``ur`` bound the link cell (from the periodic
        box when given, else the coordinate extents)."""
        if ll is None or ur is None:
            ll = coords.min(axis=0) - self.cutoff
            ur = coords.max(axis=0) + self.cutoff
        LC = Linkcell(np.array([ll, ur]), self.cutoff)
        # bucket every bond into its link cell by midpoint, vectorized.  The scalar cell
        # index uses exactly LC.ldx_of_cellndx's formula so the per-ring neighbor search
        # (LC.ldx_searchlist_of_ldx) sees the same cell numbering.
        bondlist_per_cell = {}
        if len(self._bonds) > 0:
            bond_mid = coords[self._bond_rows].mean(axis=1)
            nc = LC.cells_per_dim
            C = np.floor((bond_mid - LC.lower_left_corner) / LC.celldim).astype(int)
            C = np.clip(C, 0, nc - 1)
            bond_ldx = C[:, 2] * nc[0] * nc[1] + C[:, 1] * nc[0] + C[:, 0]
            for bi, ldx in enumerate(bond_ldx.tolist()):
                bondlist_per_cell.setdefault(ldx, []).append(bi)
        piercespecs = []
        rdict = {}
        for ring in rings_to_check:
            ring.P = coords[ring._rows]
            ring.calculate_stuff()
            oc = LC.ldx_of_cellndx(LC.cellndx_of_point(ring.COM))
            ring_atoms = set(ring.idx_list)
            for sc in LC.ldx_searchlist_of_ldx(oc):
                for bi in bondlist_per_cell.get(sc, ()):
                    bond = self._bonds[bi]
                    if ring_atoms.intersection(bond.idx_list):
                        continue
                    bond.P = coords[self._bond_rows[bi]]
                    bond.calculate_stuff()
                    test_bond = bond.mic_shift(ring.COM, box) if box is not None else bond
                    pdict = ring.pierced_by(test_bond)
                    if pdict['pierced']:
                        piercespecs.append(self._piercespec(ring, test_bond, bi))
                    if 'reason' in pdict:
                        rdict[pdict['reason']] = rdict.get(pdict['reason'], 0) + 1
        for k, v in rdict.items():
            if v:
                logger.debug(f'{k}: {v}')
        return piercespecs

    def load_coords(self, pdb):
        """Read a PDB and return its coordinates as an (Natoms x 3) array in PSF atom order
        (so ``coords[row]`` lines up with :attr:`_row_of_serial`)."""
        coorddf = coorddf_from_pdb(pdb)
        assert coorddf.shape[0] == len(self.topol.atoms), f'{pdb} is incongruent with the PSF'
        return coorddf[['x', 'y', 'z']].values

    def check_coords(self, coords, box, only_piercees):
        """Targeted pierced-ring check on an in-memory coordinate array (see
        :meth:`_check_fast`); no PDB round-trip, so a trial rotation can be scored directly."""
        return self._check_fast(coords, box, only_piercees)

    def clash_count(self, coords, moved_rows, cutoff=2.0, box=None):
        """Number of heavy-atom contacts closer than ``cutoff`` between the moved atoms
        (``moved_rows``) and the rest of the system, ignoring pairs that are covalently
        bonded across the boundary.  A cheap proxy for how sterically clean a trial pose is
        (lower is better).

        **Pass ``box``.**  Without it the scan is non-periodic, and a rotation that swings a
        tail into a neighbor ACROSS a cell face scores zero clashes while putting two heavy
        atoms on top of each other.  Measured 2026-10-16 on an example-17 patchB: a rotation
        reported as "0 new heavy-atom clash(es)" left ``SOPE 18 C29`` and ``SOPE 98 C26``
        **0.451 A** apart under the minimum image -- 65.7 A apart in raw coordinates, which is
        what this saw.  The resolver then handed that structure to dynamics.

        This is the same blindness ``_check_fast`` had before ``d89cc518``: that fix went to
        the piercing detector and not to the clash scorer beside it, which is why a module is
        not a unit of correctness here."""
        moved = np.asarray(sorted(set(int(r) for r in moved_rows)), dtype=int)
        if moved.size == 0:
            return 0
        moved_set = set(moved.tolist())
        # unmoved atoms bonded to a moved atom are expected close contacts (the hinge) -> skip
        bonded = set()
        for r in moved:
            for nb in self.topol.G.neighbors(self._serials[r]):
                nr = self._row_of_serial.get(nb)
                if nr is not None and nr not in moved_set:
                    bonded.add(nr)
        mv = moved[self._heavy[moved]]
        um = np.array([r for r in range(len(self._serials))
                       if r not in moved_set and r not in bonded and self._heavy[r]], dtype=int)
        if mv.size == 0 or um.size == 0:
            return 0
        if box is not None:
            # `box` is the 3x3 cell, as `_check_fast` and `check_coords` take it; accept a bare
            # 3-vector of edges too, so a caller holding one does not have to fake a matrix.
            L = np.asarray(box, dtype=float)
            L = np.diagonal(L).astype(float) if L.ndim == 2 else L.astype(float)
            tree = cKDTree(np.mod(coords[um], L), boxsize=L)
            probe = np.mod(coords[mv], L)
        else:
            tree = cKDTree(coords[um])
            probe = coords[mv]
        return int(sum(len(hits) for hits in tree.query_ball_point(probe, cutoff)))

    def pendant_axes(self, bond_serials, segname, max_pendant=180, max_axes=10):
        """Candidate hinge axes for swinging a glycan pendant (that contains the piercing
        bond) out of a ring.

        A ring is threaded by the bond ``bond_serials`` of glycan segment ``segname``.  To
        un-thread it, the whole sub-branch carrying that bond is rotated as a rigid body
        about an upstream rotatable bond -- a *bridge* of the glycan graph (ring bonds are
        not bridges, so they are never chosen).  For each bridge, the branch containing the
        piercing bond is the pendant.  Branches larger than ``max_pendant`` atoms are skipped
        (too much collateral motion), and the axes are returned smallest-branch-first (least
        disturbance) as dicts with the hinge endpoints' rows (``i_row`` on the axis, ``j_row``
        the pivot end of the moving branch) and ``prows`` (the moving branch's atom rows).
        """
        a, b = bond_serials
        seg_atoms = [at.serial for at in self._atoms if at.segname == segname]
        Gg = self.topol.G.subgraph(seg_atoms).copy()
        if a not in Gg or b not in Gg:
            return []
        # Restrict to the MOLECULE carrying the piercing bond, not the whole segment.  For a
        # glycan the two are the same thing -- one glycan per segment -- which is the case this
        # was written for.  A membrane puts every lipid in one segment: measured on an ex17
        # build, segment L1 holds 56,535 atoms in 538 molecules with 52,173 bridges, and the
        # loop below copies the graph once per bridge.  Scoped to the segment that is ~3 hours
        # per piercing; scoped to the molecule it is milliseconds, and the answer is identical
        # because a hinge in another molecule can never move this piercing bond.
        Gg = Gg.subgraph(nx.node_connected_component(Gg, a)).copy()
        N = Gg.number_of_nodes()
        axes = []
        for u, v in nx.bridges(Gg):
            if {u, v} == {a, b}:
                continue  # the piercing bond itself is a useless hinge
            H = Gg.copy()
            H.remove_edge(u, v)
            cu = nx.node_connected_component(H, u)
            # a and b are directly bonded, so they always share a side; the endpoint on
            # that side (inner) is the pendant we rotate, hinging on the far endpoint (outer)
            if a in cu:
                inner, outer, branch = u, v, cu
            else:
                inner, outer, branch = v, u, nx.node_connected_component(H, v)
            if len(branch) > max_pendant:
                continue
            prows = np.array([self._row_of_serial[s] for s in branch], dtype=int)
            axes.append((len(branch), dict(i_ser=outer, j_ser=inner,
                                           i_row=self._row_of_serial[outer],
                                           j_row=self._row_of_serial[inner],
                                           prows=prows)))
        axes.sort(key=lambda t: t[0])
        return [ax for _, ax in axes[:max_axes]]

    def _piercespec(self, ring, bond, bond_row):
        seg, resid = self._bond_meta[bond_row]
        return dict(
            piercer=dict(segname=seg, resid=resid,
                         resname=self.resname_of.get((seg, str(resid)), '?'),
                         segtype=self.segname_to_segtype.get(seg, 'unknown'),
                         bond_serials=list(bond.idx_list)),
            piercee=dict(segname=ring.segname, resid=ring.resid,
                         resname=self.resname_of.get((ring.segname, str(ring.resid)), '?'),
                         segtype=self.segname_to_segtype.get(ring.segname, 'unknown'),
                         ring_serials=list(ring.idx_list)),
        )

    def _check_fast(self, coords, box, only_piercees):
        """Vectorized targeted check: test only the named ``only_piercees`` rings, finding
        candidate bonds by a direct distance filter instead of a whole-system link cell.

        Exact, not approximate -- :meth:`PSFRing.pierced_by` itself rejects any bond whose
        midpoint is beyond ``cutoff`` of the ring COM, which is the same set this filter
        selects.  It avoids the per-element pandas indexing and cell assignment that make a
        full :meth:`check` cost tens of seconds.

        **The candidate filter applies the minimum image.**  It did not, and the claim of
        exactness above was false for exactly the systems this is used on.  The distance was
        taken on raw coordinates and ``mic_shift`` applied only to bonds that had already
        survived the filter, so a bond lying across a periodic boundary from the ring was
        discarded before it could be shifted.  Measured on an ex17 calibration patch
        (2026-10-13, box 76 x 66 x 86 A): a real piercing of ``CHL1 L1-199`` sat 1.13 A from
        the ring under the minimum image and **66.56 A** in raw coordinates, so the filter
        dropped it and this method reported the ring unpierced -- while the full
        :meth:`check`, which is image-aware, reported it.

        That mattered far beyond a missed report.  The rotation resolver's acceptance test is
        ``check_coords(candidate, box, [target])``; for any boundary-crossing piercing it
        returned "not pierced" whatever the candidate did, so the FIRST rotation tried was
        accepted and logged as a success while the piercing survived to be deleted later.  In
        a membrane, where lipids straddle the cell edge by construction, that is the common
        case and not an edge case.
        """
        piercespecs = []
        if len(self._bonds) == 0:
            return piercespecs
        want = {(str(s), str(r)) for s, r in only_piercees}
        target_rings = [ring for key in want for ring in self._ring_index.get(key, [])]
        if not target_rings:
            return piercespecs
        bond_mid = coords[self._bond_rows].mean(axis=1)  # (Nbonds, 3)
        # Orthorhombic cell edges for the minimum image, or None in a non-periodic run.
        edges = np.diagonal(box).astype(float) if box is not None else None
        for ring in target_rings:
            ring.P = coords[ring._rows]
            ring.calculate_stuff()
            sep = bond_mid - ring.COM
            if edges is not None:
                # Minimum image BEFORE the cutoff test, not after.  Without this a bond one
                # angstrom from the ring across the boundary measures a box-length away and is
                # filtered out, and the ring is reported unpierced.
                sep = sep - edges * np.round(sep / edges)
            d = np.linalg.norm(sep, axis=1)
            ring_atoms = set(ring.idx_list)
            for bi in np.where(d <= self.cutoff)[0]:
                bond = self._bonds[bi]
                if ring_atoms.intersection(bond.idx_list):
                    continue
                bond.P = coords[self._bond_rows[bi]]
                bond.calculate_stuff()
                test_bond = bond.mic_shift(ring.COM, box) if box is not None else bond
                if ring.pierced_by(test_bond)['pierced']:
                    piercespecs.append(self._piercespec(ring, test_bond, bi))
        return piercespecs


@countTime
def ring_check(psf, pdb, xsc=None, cutoff=4.0, segtypes=['lipid'], max_ring_size=7, only_piercees=None):
    """Convenience wrapper: build a :class:`RingChecker` and check one coordinate frame.

    For repeated checks of the *same* topology (e.g. trial side-chain rotamers), build one
    :class:`RingChecker` and call :meth:`RingChecker.check` per frame to avoid re-parsing
    the PSF and re-finding the rings each time.
    """
    return RingChecker(psf, cutoff=cutoff, segtypes=segtypes,
                       max_ring_size=max_ring_size).check(pdb, xsc=xsc, only_piercees=only_piercees)