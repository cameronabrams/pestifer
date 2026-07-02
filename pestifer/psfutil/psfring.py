#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the Ring class and the ring_check function for detection of pierced rings
"""

import logging
import networkx as nx
import numpy as np
import time

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
    
@countTime
def ring_check(psf,pdb,xsc=None,cutoff=10.0,segtypes=['lipid'],max_ring_size=7,only_piercees=None):
    """
    Check for rings in a molecular structure and identify bonds that pierce them.
    This function reads a PSF file and a PDB file, extracts the coordinates of atoms, and checks for rings in the structure.
    It identifies bonds that pierce the rings based on a specified cutoff distance and returns a list of specifications for the piercing bonds.

    Parameters
    ----------
    psf : str
        Path to the PSF file.
    pdb : str
        Path to the PDB file.
    xsc : str or None
        Path to the XSC file.  If ``None``, the system is treated as non-periodic
        (vacuum): a bounding box is derived from the coordinate extents and the
        minimum-image convention is not applied.
    cutoff : float
        Cutoff distance in Angstroms for identifying piercing bonds. Default is 10.0.
    segtypes : list
        List of segment types to include in the analysis. Default is ['lipid'].
    """
    coorddf = coorddf_from_pdb(pdb)
    if xsc is not None:
        box, orig = cell_from_xsc(xsc)
        sidelengths = np.diagonal(box)
        ll = orig - 0.5 * sidelengths
        ur = orig + 0.5 * sidelengths
    else:
        box = None
        coords = coorddf[['x', 'y', 'z']].values
        ll = coords.min(axis=0) - cutoff
        ur = coords.max(axis=0) + cutoff
        logger.debug('No XSC file — treating system as non-periodic (vacuum)')
    LC = Linkcell(np.array([ll, ur]), cutoff)
    topol = PSFContents(psf, parse_topology=['bonds'], topology_segtypes=segtypes)
    assert coorddf.shape[0] == len(topol.atoms), f'{psf} and {pdb} are incongruent'
    coorddf['segname'] = [a.segname for a in topol.atoms.data]
    segname_to_segtype = {a.segname: a.segtype for a in topol.atoms}
    # key on str(resid): PSF atoms carry an int resid but the ring/bond resid comes from
    # the ingested coordinate frame as a string, so the tuples must be normalized to match
    resname_of = {(a.segname, str(a.resid)): a.resname for a in topol.atoms}
    logger.debug(f'ingesting coords into bonds...(could take a while)')
    topol.bonds.ingest_coordinates(coorddf, pos_key=['x', 'y', 'z'], meta_key=['segname', 'resid'])
    logger.debug(f'Assiging link-cell indices to each bond')
    topol.bonds.assign_cell_indexes(LC)
    logger.debug(f'Sorting bonds by link-cell index')
    bondlist_per_cell = {}
    for b in topol.bonds:
        if not b.linkcell_idx in bondlist_per_cell:
            bondlist_per_cell[b.linkcell_idx] = PSFBondList([])
        bondlist_per_cell[b.linkcell_idx].append(b)
    Rings = RingList(topol.G, length_bound=max_ring_size)
    logger.debug(f'{len(Rings)} rings to be analyzed for piercings')
    logger.debug(f'ingesting coords into rings...')
    Rings.ingest_coordinates(coorddf, pos_key=['x', 'y', 'z'], meta_key=['segname', 'resid'])
    if only_piercees is not None:
        # targeted re-check: keep only the named (segname, resid) rings, so re-checking
        # one candidate side-chain rotation does not re-scan every ring in the system
        want = {(str(s), str(r)) for s, r in only_piercees}
        Rings.data = [ring for ring in Rings.data if (str(ring.segname), str(ring.resid)) in want]
        logger.debug(f'restricting to {len(Rings)} ring(s) matching only_piercees')
    piercespecs = []
    rdict = {}
    dtsum = 0.0
    nr = 0
    logger.debug(f'Examining {len(Rings)} rings...')
    for iring, ring in enumerate(Rings):
        t1 = time.time()
        # logger.debug(f'Ring {iring} {ring.segname} {ring.resid}...')
        oc = LC.ldx_of_cellndx(LC.cellndx_of_point(ring.COM))
        # logger.debug(f'ring {str(ring)} center in cell {oc} -- members in {ring.M["linkcell_idx"]}')
        searchcells = LC.ldx_searchlist_of_ldx(oc)
        search_bonds = PSFBondList([])
        # logger.debug(f'searchcells {searchcells}')
        for sc in searchcells:
            search_bonds.extend(bondlist_per_cell.get(sc, PSFBondList([])))
        # logger.debug(f'Ring {iring+1}/{len(Rings)} in c {oc} searching c {searchcells} captures {len(search_bonds)} bonds')
        piercing_bonds = PSFBondList([])
        for bond in search_bonds:
            if any([x in ring.idx_list for x in bond.idx_list]):
                continue
            test_bond = bond.mic_shift(ring.COM, box) if box is not None else bond
            pdict = ring.pierced_by(test_bond)
            if pdict['pierced']: piercing_bonds.append(test_bond)
            if 'reason' in pdict:
                if not pdict['reason'] in rdict:
                    rdict[pdict['reason']] = 0
                rdict[pdict['reason']] += 1
                if 'error' in pdict:
                    if not 'error' in rdict:
                        rdict['error'] = []
                    if pdict['error'] < 1.0:
                        rdict['error'].append(dict(pdict=pdict, bond=bond.idx_list, ring=ring.idx_list))
        t2 = time.time()
        dt = t2 - t1
        dtsum += dt
        nr += 1
        dtav = dtsum / nr
        if len(Rings) < 10 or iring % 10 == 0: logger.debug(f'{dtav:.3f} s per ring, {nr}/{len(Rings)} rings, {dtav * (len(Rings) - nr):.3f} s remaining...')
        if len(piercing_bonds) > 0:
            ringname = '-|- ' + ' -- '.join([str(topol.included_atoms.get((lambda a: a.serial == x))) for x in ring.idx_list]) + ' -|-'
            logger.debug(f'ring {ringname}:')
            for bond in piercing_bonds:
                piercespecs.append(dict(
                    piercer=dict(segname=bond.segname, resid=bond.resid, resname=resname_of.get((bond.segname, str(bond.resid)), '?'), segtype=segname_to_segtype.get(bond.segname, 'unknown')),
                    piercee=dict(segname=ring.segname, resid=ring.resid, resname=resname_of.get((ring.segname, str(ring.resid)), '?'), segtype=segname_to_segtype.get(ring.segname, 'unknown')),
                ))
                bondname = ' -- '.join([str(topol.included_atoms.get((lambda a: a.serial == x))) for x in bond.idx_list])
                logger.debug(f'  pierced by bond [ {bondname} ]')
    for k, v in rdict.items():
        if v:
            if type(v) == list:
                for e in v:
                    if type(e) == dict:
                       for kk, vv in e.items():
                            logger.debug(f'    {kk}: {vv}')
                    else:
                       logger.debug(f'{k}: {v}')
            else:
                logger.debug(f'{k}: {v}')
    return piercespecs