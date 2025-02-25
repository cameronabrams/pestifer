#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the Ring class and the ring_check function for detection of pierced rings
"""

import logging
import networkx as nx
import numpy as np
import time

from functools import singledispatchmethod
from itertools import pairwise

from .psfutil.psfbond import PSFBondList
from .psfutil.psfcontents import PSFContents
from .psfutil.psftopoelement import PSFTopoElementList,PSFTopoElement

from .util.coord import coorddf_from_pdb
from .util.linkcell import Linkcell
from .util.util import countTime, cell_from_xsc

logger=logging.getLogger(__name__)

def lawofcos(a,b):
    """lawofcos return the cosine of the angle defined by vectors a and b if they share a vertex (the LAW OF COSINES)

    :param a: a vector
    :type a: numpy.ndarray(3,float)
    :param b: another vector
    :type b: numpy.ndarray(3,float)
    :return: cosine of the angle formed by a and b
    :rtype: float
    """
    return np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b))

class Ring(PSFTopoElement):

    def treadmill(self):
        """ yield the treadmilled versions of the list """
        for i in range(1,len(self.idx_list)):
            yield self.idx_list[i:]+self.idx_list[:i]

    def __eq__(self,other):
        check1=any([self.idx_list==other.idx_list] + [self.idx_list==x for x in other.treadmill()])
        check2=any([self.idx_list==other.idx_list[::-1]]+[self.idx_list==x[::-1] for x in other.treadmill()])
        return check1 or check2

    # @countTime
    def pierced_by(self,bond,cutoff=3.5,tol=1.e-5):
        """determines if PSFBond B pierces ring
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
    @singledispatchmethod
    def __init__(self,input_obj):
        self.data=input_obj
    @__init__.register(nx.Graph)
    def _from_graph(self,G):
        L=[]
        for ll in nx.chordless_cycles(G):
            L.append(Ring(ll))
        super().__init__(L)

    def __str__(self):
        return ';'.join([str(x) for x in self])
    
@countTime
def ring_check(psf,pdb,xsc,cutoff=10.0,segtypes=['lipid']):
    box,orig=cell_from_xsc(xsc)
    coorddf=coorddf_from_pdb(pdb)
    sidelengths=np.diagonal(box)
    ll=orig-0.5*sidelengths
    ur=orig+0.5*sidelengths
    LC=Linkcell(np.array([ll,ur]),cutoff)
    topol=PSFContents(psf,parse_topology=['bonds'],topology_segtypes=segtypes)
    assert coorddf.shape[0]==len(topol.atoms),f'{psf} and {pdb} are incongruent'
    coorddf['segname']=[a.chainID for a in topol.atoms]
    logger.debug(f'Injesting coords into bonds...(could take a while)')
    topol.bonds.injest_coordinates(coorddf,pos_key=['x','y','z'],meta_key=['segname','resid'])
    logger.debug(f'Assiging link-cell indices to each bond')
    topol.bonds.assign_cell_indexes(LC)
    logger.debug(f'Sorting bonds by link-cell index')
    bondlist_per_cell={}
    for b in topol.bonds:
        if not b.linkcell_idx in bondlist_per_cell:
            bondlist_per_cell[b.linkcell_idx]=PSFBondList([])
        bondlist_per_cell[b.linkcell_idx].append(b)
    Rings=RingList(topol.G)
    logger.debug(f'{len(Rings)} rings to be analyzed for piercings')
    logger.debug(f'Injesting coords into rings...')
    Rings.injest_coordinates(coorddf,pos_key=['x','y','z'],meta_key=['segname','resid'])
    piercespecs=[]
    rdict={}
    dtsum=0.0
    nr=0
    logger.debug(f'Examining {len(Rings)} rings...')
    for iring,ring in enumerate(Rings):
        t1=time.time()
        # logger.debug(f'Ring {iring} {ring.segname} {ring.resid}...')
        oc=LC.ldx_of_cellndx(LC.cellndx_of_point(ring.COM))
        # logger.debug(f'ring {str(ring)} center in cell {oc} -- members in {ring.M["linkcell_idx"]}')
        searchcells=LC.ldx_searchlist_of_ldx(oc)
        search_bonds=PSFBondList([])
        # logger.debug(f'searchcells {searchcells}')
        for sc in searchcells:
            search_bonds.extend(bondlist_per_cell.get(sc,PSFBondList([])))
        # logger.debug(f'Ring {iring+1}/{len(Rings)} in c {oc} searching c {searchcells} captures {len(search_bonds)} bonds')
        piercing_bonds=PSFBondList([])
        for bond in search_bonds:
            if any([x in ring.idx_list for x in bond.idx_list]):
                continue
            test_bond=bond.mic_shift(ring.COM,box)
            pdict=ring.pierced_by(test_bond)
            if pdict['pierced']: piercing_bonds.append(test_bond)
            if 'reason' in pdict:
                if not pdict['reason'] in rdict:
                    rdict[pdict['reason']]=0
                rdict[pdict['reason']]+=1
                if 'error' in pdict:
                    if not 'error' in rdict:
                        rdict['error']=[]
                    if pdict['error']<1.0:
                        rdict['error'].append(dict(pdict=pdict,bond=bond.idx_list,ring=ring.idx_list))
        t2=time.time()
        dt=t2-t1
        dtsum+=dt
        nr+=1
        dtav=dtsum/nr
        if len(Rings)<10 or iring%10==0: logger.debug(f'{dtav:.3f} s per ring, {nr}/{len(Rings)} rings, {dtav*(len(Rings)-nr):.3f} s remaining...')
        if len(piercing_bonds)>0:
            ringname='-|- '+ ' -- '.join([str(topol.included_atoms.get(serial=x)[0]) for x in ring.idx_list]) + ' -|-'
            logger.debug(f'ring {ringname}:')
            for bond in piercing_bonds:
                piercespecs.append(dict(piercer=dict(segname=bond.segname,resid=bond.resid),piercee=dict(segname=ring.segname,resid=ring.resid)))
                bondname=' -- '.join([str(topol.included_atoms.get(serial=x)[0]) for x in bond.idx_list])
                logger.debug(f'  pierced by bond [ {bondname} ]')
    for k,v in rdict.items():
        if v:
            if type(v)==list:
                for e in v:
                    if type(e)==dict:
                       for kk,vv in e.items():
                            logger.debug(f'    {kk}: {vv}')
                    else:
                       logger.debug(f'{k}: {v}')
            else:
                logger.debug(f'{k}: {v}')
    return piercespecs