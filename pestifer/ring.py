#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Rings -- checks for pierced rings based on topology and coordinate information
"""

import numpy as np
import os
from collections import UserList
from functools import singledispatchmethod
import networkx as nx
from copy import deepcopy
import logging
logger=logging.getLogger(__name__)
from .psf import PSFContents,PSFBondList,PSFTopoElementList,PSFTopoElement
from .linkcell import Linkcell
from .util import cell_from_xsc
from .coord import coorddf_from_pdb, mic_shift
from itertools import compress

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

    def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ']):
        super().injest_coordinates(A,idx_key=idx_key,pos_key=pos_key)
        # logger.debug(f'COM {self.COM}')
        iR=Ring(list(range(len(self.idx_list))))
        a=iR.treadmill()
        self.B=[]
        for i,j in zip(iR.idx_list,next(a)):
            b=self.P[i]-self.P[j]
            # logger.debug(f'R-{self.idx[i]}-{self.idx[j]}: {b}')
            self.B.append(b)
            # logger.debug(f'R-{self.idx[i]}-{self.idx[j]}: {self.B[-1]}')
        self.B=np.array(self.B)
        # logger.debug(f'B {self.B}')
        # Approximate the normal unit vector of the plane of the ring
        # by averaging cross-products of all adjacent origin-to-vertex
        # vectors.  This orients the unit vector such that it is 
        # in the +z direction in the local coordinate system defined
        # by the ring with vertices ordered in the counter-clockwise
        # direction in the plane.
        iR=Ring(list(range(len(self.idx_list))))
        a=iR.treadmill()
        self.C=[]
        for i,j in zip(iR.idx_list,next(a)):
            c=np.cross(self.B[i],self.B[j])
            # logger.debug(f'C-{self.idx[i]}-{self.idx[j]}: {c}')
            self.C.append(c)
        self.C=np.array(self.C)
        # logger.debug(f'C {self.C}')
        n=np.sum(self.C,axis=0)
        self.n=n/np.linalg.norm(n)
        # logger.debug(f'n {self.n}')
        # compute planarity as average of all
        # cross-i-cross-i+1 dot products
        iR=Ring(list(range(len(self.idx_list))))
        a=iR.treadmill()
        self.planarity=0
        for i,j in zip(iR.idx_list,next(a)):
            ci=self.C[i]
            cj=self.C[j]
            self.planarity+=lawofcos(ci,cj)
        self.planarity/=len(self)
        # get the d for n-dot-r + d = 0 equation of the plane
        # n[0]*(x-O[0])+n[1]*(y-O[1])+n[2]*(z-O[2])=0
        # n[0]*x + n[1]*y + n[2]*z - (n[0]*O[0]+n[1]*O[1]+n[2]*O[2]) = 0
        self.d=-np.dot(self.n,self.COM)        

    def treadmill(self):
        """ yield the treadmilled versions of the list """
        for i in range(1,len(self.idx_list)):
            yield self.idx_list[i:]+self.idx_list[:i]

    def __eq__(self,other):
        check1=any([self.idx_list==other.idx_list] + [self.idx_list==x for x in other.treadmill()])
        check2=any([self.idx_list==other.idx_list[::-1]]+[self.idx_list==x[::-1] for x in other.treadmill()])
        return check1 or check2

    def pierced_by(self,B,tol=1.e-2):
        """determines if PSFBond B pierces ring; 
        uses ray projection method and fact that scaled length must be between 
        0 and 1 for a plane intersection
        """
        V=B.P[0]-B.P[1]
        # If these two points are on either side of the ring plane,
        # then the length of the normal vector from P[0] to the plane
        # will be less than the length of the vector from P[1] to P[0] projected
        # onto the plane normal. 
        denom=np.linalg.norm(np.dot(self.n,V)*self.n) # length of vector from P[1] to P[0] projected onto normal
        OP=np.dot(self.n,B.P[0]-self.COM)*self.n # normal vector from plane to P[0]
        numer=np.linalg.norm(OP) # length of normal vector from plane to P[0]
        # logger.debug(f'numer {numer} denom {denom}')
        t=-1.0
        if denom>1.e-13:
            t=numer/denom # a fraction if P[0] and P[1] are on either side of the plane
        # logger.debug(f'{str(self)} b {P} t {t}')
        if 0<t<1:
            # compute point in ring plane that marks intersection with this vector
            # V=P[0]-P[1]
            # P[1]=P[0]-V
            # intersection point = P[0]-t*V
            PP=B.P[0]-t*V
            # logger.debug(f't {t} PP {PP}')
            # determine if PP is inside ring:
            # compute the series of unit-vector cross-products v(i)-PP-v((i+1)%N)
            # sum will have a large component along normal vector if yes, essentially 0 if no
            iR=Ring(list(range(len(self.idx_list))))
            a=iR.treadmill()
            sumC=np.zeros(3)
            for i,j in zip(iR.idx_list,next(a)):
                V1=PP-self.P[i]
                V2=PP-self.P[j]
                c=np.cross(V1/np.linalg.norm(V1),V2/np.linalg.norm(V2))
                sumC+=c/np.linalg.norm(c)
            tst=np.dot(self.n,sumC/len(self))
            is_inside=np.abs(tst-1.0)<tol
            # logger.debug(f'tst {tst} is_inside {is_inside}')
            return is_inside,PP
        return False,np.ones(3)*np.nan

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
    

def ring_check(psf,pdb,xsc):
    box,orig=cell_from_xsc(xsc)
    coorddf=coorddf_from_pdb(pdb)
    sidelengths=np.diagonal(box)
    ll=orig-0.5*sidelengths
    ur=orig+0.5*sidelengths
    LC=Linkcell(np.array([ll,ur]),10.0,atidxlabel=None)
    LC.populate(coorddf,os.cpu_count())
    topol=PSFContents(psf,parse_topology=['bonds'])
    topol.nonsolvent_atoms.injest_coordinates(coorddf,idx_key=None,pos_key=['x','y','z'])
    topol.bonds.injest_coordinates(coorddf,idx_key=None,pos_key=['x','y','z'])
    topol.bonds.validate_images(box)
    for bond in topol.bonds:
        bond.oc_ldx=LC.ldx_of_cellndx(LC.cellndx_of_point(bond.COM))
    Rings=RingList(topol.G)
    logger.debug(f'{len(Rings)} rings to be analyzed for piercings')
    Rings.injest_coordinates(coorddf,idx_key=None,pos_key=['x','y','z'])
    piercing_tally=[]
    piercespecs=[]
    for iring,ring in enumerate(Rings):
        iat=topol.nonsolvent_atoms.get(serial=ring.idx_list[0])[0]
        chain,resid=iat.chainID,iat.resseqnum
        oc=LC.ldx_of_cellndx(LC.cellndx_of_point(ring.COM))
        # logger.debug(f'ring {str(ring)} center in cell {oc} -- members in {ring.M["linkcell_idx"]}')
        nc=LC.searchlist_of_ldx(oc)
        # logger.debug(f'...search for bonds in cells {nc} and {oc}')
        search_bonds=PSFBondList([])
        for ldx in nc+[oc]:
            this_sb=[b for b in topol.bonds if b.oc_ldx==ldx]
            search_bonds.extend(this_sb)
        # logger.debug(f'...search will include {len(search_bonds)} bonds')
        piercings=[]
        for bond in search_bonds:
            bc=bond.copy().mic_shift(ring.COM,box)
            p,r=ring.pierced_by(bc)
            piercings.append(p)
        piercing_bonds=list(compress(search_bonds,piercings))
        if len(piercing_bonds)>0:
            ringname='- '+ ' -- '.join([str(topol.nonsolvent_atoms.get(serial=x)) for x in ring.idx_list]) + ' -'
            logger.debug(f'ring {ringname}:')
            for b in piercing_bonds:
                ai=topol.nonsolvent_atoms.get(serial=b.serial1)[0]
                aj=topol.nonsolvent_atoms.get(serial=b.serial2)[0]
                pchain,presid=ai.chainID,ai.resseqnum
                if chain==pchain and resid==presid:
                    continue
                # logger.debug(f'...{b.serial1}--{b.serial2}')
                piercespecs.append(dict(piercee=dict(chain=chain,resid=resid),piercer=dict(chain=pchain,resid=presid)))
                logger.debug(f'  pierced by bond [ {str(ai)} -- {str(aj)} ]')
        piercing_tally.append(any(piercings))
        logger.debug(f'ring {iring}')
    return piercespecs