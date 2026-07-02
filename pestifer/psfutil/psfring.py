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
        - The COM of the bond is not more than half the distance between the two endpoints of the bond.
        - The ring's points, when projected onto a plane perpendicular to the bond at the bond's midpoint, form a closed loop (i.e., the sum of angles around the COM is approximately 2π).

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
        # - reject if bond COM and ring COM are more than the cutoff from each other
        # - reject if both bond endpoints are on the same side of the plane
        dist=np.linalg.norm(bond.COM-self.COM)
        # logger.debug(f'{self.COM} {bond.COM} {dist:.3f} {cutoff} {0.5*bond.b[0,1]}')
        if dist>cutoff:
            return dict(pierced=False,reason='cutoff')
        elif dist>0.5*bond.b[0,1]:
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
        Cutoff distance in Angstroms for identifying piercing bonds.
    segtypes : list
        Segment types whose rings are checked.
    max_ring_size : int
        Only rings of at most this many atoms are considered.
    """

    def __init__(self, psf, cutoff=10.0, segtypes=['lipid'], max_ring_size=7):
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
        constructor are reused.  ``xsc`` sets the periodic box (``None`` -> vacuum, box from
        coordinate extents).  ``only_piercees`` (iterable of ``(segname, resid)``) restricts
        which rings are tested -- a cheap targeted re-check after a trial rotation.
        """
        topol = self.topol
        coorddf = coorddf_from_pdb(pdb)
        assert coorddf.shape[0] == len(topol.atoms), f'{pdb} is incongruent with the PSF'
        box = cell_from_xsc(xsc)[0] if xsc is not None else None
        if only_piercees is not None:
            # targeted re-check after a trial rotation: vectorized, no link cell.  Only a
            # handful of named rings are tested, so skip the whole-system coordinate ingest.
            coords = coorddf[['x', 'y', 'z']].values
            return self._check_fast(coords, box, only_piercees)
        if xsc is not None:
            orig = cell_from_xsc(xsc)[1]
            sidelengths = np.diagonal(box)
            ll = orig - 0.5 * sidelengths
            ur = orig + 0.5 * sidelengths
        else:
            ll = None
            logger.debug('No XSC file — treating system as non-periodic (vacuum)')
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
            bond_ldx = C[:, 2] * nc[0] * nc[1] + C[:, 1] * nc[1] + C[:, 0]
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

    def clash_count(self, coords, moved_rows, cutoff=2.0):
        """Number of heavy-atom contacts closer than ``cutoff`` between the moved atoms
        (``moved_rows``) and the rest of the system, ignoring pairs that are covalently
        bonded across the boundary.  A cheap proxy for how sterically clean a trial pose is
        (lower is better)."""
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
        tree = cKDTree(coords[um])
        return int(sum(len(hits) for hits in tree.query_ball_point(coords[mv], cutoff)))

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

        This is exact, not approximate -- :meth:`PSFRing.pierced_by` itself rejects any bond
        whose midpoint is beyond ``cutoff`` of the ring COM, which is the same set the direct
        filter selects.  It avoids the per-element pandas indexing and cell assignment that
        make a full :meth:`check` cost tens of seconds.
        """
        piercespecs = []
        if len(self._bonds) == 0:
            return piercespecs
        want = {(str(s), str(r)) for s, r in only_piercees}
        target_rings = [ring for key in want for ring in self._ring_index.get(key, [])]
        if not target_rings:
            return piercespecs
        bond_mid = coords[self._bond_rows].mean(axis=1)  # (Nbonds, 3)
        for ring in target_rings:
            ring.P = coords[ring._rows]
            ring.calculate_stuff()
            d = np.linalg.norm(bond_mid - ring.COM, axis=1)
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
def ring_check(psf, pdb, xsc=None, cutoff=10.0, segtypes=['lipid'], max_ring_size=7, only_piercees=None):
    """Convenience wrapper: build a :class:`RingChecker` and check one coordinate frame.

    For repeated checks of the *same* topology (e.g. trial side-chain rotamers), build one
    :class:`RingChecker` and call :meth:`RingChecker.check` per frame to avoid re-parsing
    the PSF and re-finding the rings each time.
    """
    return RingChecker(psf, cutoff=cutoff, segtypes=segtypes,
                       max_ring_size=max_ring_size).check(pdb, xsc=xsc, only_piercees=only_piercees)