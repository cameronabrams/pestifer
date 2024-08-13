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
from .util import cell_from_xsc,countTime
from .coord import coorddf_from_pdb
import time

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

    # def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ']):
    #     super().injest_coordinates(A,idx_key=idx_key,pos_key=pos_key)
        # logger.debug(f'COM {self.COM}')
        # iR=Ring(list(range(len(self.idx_list))))
        # a=iR.treadmill()
        # self.B=[]
        # for i,j in zip(iR.idx_list,next(a)):
        #     b=self.P[i]-self.P[j]
        #     # logger.debug(f'R-{self.idx[i]}-{self.idx[j]}: {b}')
        #     self.B.append(b)
        #     # logger.debug(f'R-{self.idx[i]}-{self.idx[j]}: {self.B[-1]}')
        # self.B=np.array(self.B)
        # # logger.debug(f'B {self.B}')
        # # Approximate the normal unit vector of the plane of the ring
        # # by averaging cross-products of all adjacent origin-to-vertex
        # # vectors.  This orients the unit vector such that it is 
        # # in the +z direction in the local coordinate system defined
        # # by the ring with vertices ordered in the counter-clockwise
        # # direction in the plane.
        # iR=Ring(list(range(len(self.idx_list))))
        # a=iR.treadmill()
        # self.C=[]
        # for i,j in zip(iR.idx_list,next(a)):
        #     c=np.cross(self.B[i],self.B[j])
        #     # logger.debug(f'C-{self.idx[i]}-{self.idx[j]}: {c}')
        #     self.C.append(c)
        # self.C=np.array(self.C)
        # # logger.debug(f'C {self.C}')
        # n=np.sum(self.C,axis=0)
        # self.n=n/np.linalg.norm(n)
        # # logger.debug(f'n {self.n}')
        # # compute planarity as average of all
        # # cross-i-cross-i+1 dot products
        # iR=Ring(list(range(len(self.idx_list))))
        # a=iR.treadmill()
        # self.planarity=0
        # for i,j in zip(iR.idx_list,next(a)):
        #     ci=self.C[i]
        #     cj=self.C[j]
        #     self.planarity+=lawofcos(ci,cj)
        # self.planarity/=len(self)
        # # get the d for n-dot-r + d = 0 equation of the plane
        # # n[0]*(x-O[0])+n[1]*(y-O[1])+n[2]*(z-O[2])=0
        # # n[0]*x + n[1]*y + n[2]*z - (n[0]*O[0]+n[1]*O[1]+n[2]*O[2]) = 0
        # self.d=-np.dot(self.n,self.COM)        

    def treadmill(self):
        """ yield the treadmilled versions of the list """
        for i in range(1,len(self.idx_list)):
            yield self.idx_list[i:]+self.idx_list[:i]

    def __eq__(self,other):
        check1=any([self.idx_list==other.idx_list] + [self.idx_list==x for x in other.treadmill()])
        check2=any([self.idx_list==other.idx_list[::-1]]+[self.idx_list==x[::-1] for x in other.treadmill()])
        return check1 or check2

    # @countTime
    def pierced_by(self,bond,cutoff=3.5,tol=1.e-3):
        """determines if PSFBond B pierces ring
        """

        # reject if bond COM and ring COM are more than the cutoff from each other
        # COM=bsrs[['comx','comy','comz']].to_numpy()
        dist=np.linalg.norm(bond.COM-self.COM)
        if dist>cutoff:
            return dict(pierced=False,reason='cutoff')

        # logger.debug(f'PIERCE: B {B.idx_list} R {self.idx_list}')
        
        # project all ring points into plane perpendicular to bond at bond midpoint
        # lbv,rbv=bsrs[['nvpx','nvpy','nvpz']].to_numpy(),bsrs[['nvmx','nvmy','nvmz']].to_numpy()
        lbv,rbv=bond.radii[0],bond.radii[1]
        # logger.debug(f'{lbv}')
        (a,b,c),(aa,bb,cc)=lbv,rbv
        d,e,f=bond.COM
        # logger.debug(f'b-com: {[d,e,f]} b-nv {[a,b,c]} {np.linalg.norm(np.array([a,b,c]))} r-com {self.COM} -- d {dist}')
        projected_ring=[]
        for p in self.P:
            x,y,z=p
            t=a*d-a*x+b*e-b*y+c*f-c*z
            # logger.debug(f'{[x,y,z]} -> {t}')
            projected_ring.append(np.array([x+t*a,y+t*b,z+t*c]))
        projected_ring=np.array(projected_ring)
        # logger.debug(f'proj ring {projected_ring}')
        pCOM=np.mean(projected_ring,axis=0)
        # reject if both bond endpoints are on the same side of the plane
        dd=np.linalg.norm(self.COM-pCOM)
        if dd>0.5*np.linalg.norm(bond.P[0]-bond.P[1]):
            return dict(pierced=False,reason='both bond atoms on same side of ring plane')

        phisum=0.0
        for i in range(len(projected_ring)):
            if i+1==len(projected_ring):
                ip1=0
            else:
                ip1=i+1
            p0=projected_ring[i]
            p1=projected_ring[ip1]
            b0=p0-bond.COM
            b1=p1-bond.COM
            # logger.debug(f'{i}:{ip1} {b0} {b1}')
            phi=np.arccos(lawofcos(b0,b1))
            phisum+=phi
            # logger.debug(f'{phi} -> {phisum}')
        diff=phisum-np.pi*2
        # logger.debug(f'{phisum} {diff} {tol}')
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
    LC=Linkcell(np.array([ll,ur]),cutoff,atidxlabel=None)
    # ncpu=os.cpu_count() if coorddf.shape[0]>1000 else 1
    # LC.populate(coorddf,ncpu,report_info=False,profile=True)
    topol=PSFContents(psf,parse_topology=['bonds'],topology_segtypes=segtypes)
    # topol.atoms.make_df()
    coorddf['segname']=[a.chainID for a in topol.atoms]#.df['chainID']
    topol.bonds.injest_coordinates(coorddf,pos_key=['x','y','z'],meta_key=['segname','resid'])
    logger.debug(f'...bond coords')
    topol.bonds.assign_cell_indexes(LC)
    logger.debug(f'...bond cell indexes')
    Rings=RingList(topol.G)
    logger.debug(f'{len(Rings)} rings to be analyzed for piercings')
    Rings.injest_coordinates(coorddf,pos_key=['x','y','z'],meta_key=['segname','resid'])
    logger.debug(f'...ring coords')
    # piercing_tally=[]
    piercespecs=[]
    bonds_in={}
    rdict={}
    dtsum=0.0
    nr=0
    for iring,ring in enumerate(Rings):
        t1=time.time()
        # logger.debug(f'Ring {iring} {ring.segname} {ring.resid}...')
        oc=LC.ldx_of_cellndx(LC.cellndx_of_point(ring.COM))
        # logger.debug(f'ring {str(ring)} center in cell {oc} -- members in {ring.M["linkcell_idx"]}')
        searchcells=LC.searchlist_of_ldx(oc)+[oc]
        # logger.debug(f'searchcells {searchcells}')
        if not oc in bonds_in:
            bonds_in[oc]=PSFBondList([])
            for sc in searchcells:
                bonds_in[oc].extend(topol.bonds.get(linkcell_idx=sc))
            # bonds_in[oc]=topol.bonddf.df[topol.bonds.df['linkcell_idx'].isin(searchcells)]
            # logger.debug(f'cell {oc} has {len(bonds_in[oc])} bonds')
        # logger.debug(f'...search for bonds in cells {nc} and {oc}')
        search_bonds=bonds_in[oc]
        # logger.debug(f'...search nominally include {len(search_bonds)} bonds')
        # search_bonds['piercer']=[False]*search_bonds.shape[0]
        piercing_bonds=PSFBondList([])
        for bond in search_bonds:
        # for i in search_bonds.index:
            # if ring.segname==bond.segname and ring.resid==bond.resid:
            #     keep[i]=False
            # elif any([x in ring.idx_list for x in bond.idx_list]):
            # if not any([x in ring.idx_list for x in search_bonds.loc[i,['i','j']]]):
            if any([x in ring.idx_list for x in bond.idx_list]):
                continue
                # idx,jdx=search_bonds.loc[i,['i','j']]
                # ri,bli=mic_shift(coorddf.loc[idx,['x','y','z']].to_numpy(),ring.COM,box)
                # ri,bli=mic_shift(bond.P[0],ring.COM,box)
                # search_bonds.loc[i,['ix','iy','iz']]=ri
                # rj,blj=mic_shift(coorddf.loc[jdx,['x','y','z']].to_numpy(),ring.COM,box)
                # rj,blj=mic_shift(bond.P[1],ring.COM,box)
                # search_bonds.loc[i,['jx','jy','jz']]=rj
                # com,blc=mic_shift(bond.COM,ring.COM,box)
                # search_bonds.loc[i,['comx','comy','comz']]=com
                # logger.debug(f'{i} {ri} {com}')
                # vp=ri-com
                # logger.debug(f'{i} {ri} {com} {vp}')
                # search_bonds.loc[i,['nvpx','nvpy','nvpz']]=vp/np.linalg.norm(vp)
                # vm=rj-com
                # search_bonds.loc[i,['nvmx','nvmy','nvmz']]=vm/np.linalg.norm(vm)
            test_bond=bond.mic_shift(ring.COM,box)
            pdict=ring.pierced_by(test_bond)
            if pdict['pierced']: piercing_bonds.append(test_bond)
            # pdict=ring.pierced_by(search_bonds.loc[i])
            # search_bonds.loc[i,'piercer']=pdict['pierced']
            if 'reason' in pdict:
                if not pdict['reason'] in rdict:
                    rdict[pdict['reason']]=0
                rdict[pdict['reason']]+=1
                if 'error' in pdict:
                    if not 'error' in rdict:
                        rdict['error']=[]
                    if pdict['error']<1.0:
                        rdict['error'].append(dict(error=pdict['error'],bond=bond.idx_list,ring=ring.idx_list))
        t2=time.time()
        dt=t2-t1
        dtsum+=dt
        nr+=1
        dtav=dtsum/nr
        if len(Rings)<10 or iring%10==0: logger.debug(f'{dtav:.3f} s per ring, {nr}/{len(Rings)} rings, {dtav*(len(Rings)-nr):.3f} s remaining...')
        # search_bonds['piercer']=piercings
        # piercing_bonds=search_bonds[search_bonds['piercer']==True]
        if len(piercing_bonds)>0:
            ringname='-|- '+ ' -- '.join([str(topol.included_atoms.get(serial=x)[0]) for x in ring.idx_list]) + ' -|-'
            logger.debug(f'ring {ringname}:')
            for bond in piercing_bonds:
                piercespecs.append(dict(piercer=dict(segname=bond.segname,resid=bond.resid),piercee=dict(segname=ring.segname,resid=ring.resid)))
                bondname=' -- '.join([str(topol.included_atoms.get(serial=x)[0]) for x in bond.idx_list])
                logger.debug(f'  pierced by bond [ {bondname} ]')
        # piercing_tally.append(any(piercings))
        # logger.debug(f'ring {iring}')
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