# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""PSF file processing
"""
from .basemod import *
from .mods import SSBond,SSBondList,Link,LinkList
from .coord import mic_shift
from .util import countTime
from .stringthings import split_ri, my_logger
from .config import segtype_of_resname, Config
import networkx as nx 
import numpy as np
import pandas as pd
import logging
from copy import deepcopy
from argparse import Namespace
from functools import singledispatchmethod
from itertools import product

logger=logging.getLogger(__name__)

class LineList(UserList):
    def __init__(self,data):
        super().__init__(data)

class PSFTopoElement(Namespace):
    def __init__(self,idx_list):
        a_dict={}
        for i,idx in enumerate(idx_list):
            a_dict.update({f'serial{i+1}':idx})
        super().__init__(**a_dict)
        self.idx_list=idx_list.copy()

    def __len__(self):
        return len(self.idx_list)
    
    def __str__(self):
        return '-'.join([str(x) for x in self.idx_list])

    def calculate_stuff(self):
        self.COM=np.mean(self.P,axis=0)
        self.Radii=self.P-self.COM
        self.radii=np.array([n/d for n,d in zip(self.Radii,np.linalg.norm(self.Radii,axis=1))])
        self.B=np.array([p[1]-p[0] for p in product(self.P,repeat=2)]).reshape((len(self.P),len(self.P),3))
        self.b=np.linalg.norm(self.B.reshape((len(self.P)**2,3)),axis=1).reshape(len(self.P),len(self.P))

    def injest_coordinates(self,A,pos_key=['posX','posY','posZ'],meta_key=[]):
        self.P=np.array(A.loc[self.idx_list][pos_key].values)
        if np.any(np.isnan(self.P)):
            logger.debug(f'nan detected in {self.P}')
        if len(self.P)>1:
            self.calculate_stuff()
        for k in meta_key:
            self.__dict__[k]=A.loc[self.idx_list[0],k]

    def assign_cell_index(self,LC):
        self.linkcell_idx=LC.ldx_of_cellndx(LC.cellndx_of_point(self.COM))

    def validate_image(self,box):
        boxlen=np.array([box[i][i] for i in [0,1,2]])
        in_same_img=np.linalg.norm(self.P-self.COM)<(boxlen/2)
        self.same_image=np.all(in_same_img)
        return self.same_image
    
    def mic_shift(self,reference_point,box):
        holder=self.copy()
        for i in range(len(self.P)):
            holder.P[i],bl=mic_shift(self.P[i],reference_point,box)
        if np.any(bl): holder.calculate_stuff()
        return holder
    
    def copy(self):
        return deepcopy(self)

    def match(self,**crits):
        return all([self.__dict__[k]==v for k,v in crits.items()])

class PSFTopoElementList(UserList):
    def __init__(self,data):
        super().__init__(data)

    @countTime
    def injest_coordinates(self,A,pos_key=['posX','posY','posZ'],meta_key=[]):
        for b in self:
            b.injest_coordinates(A,pos_key=pos_key,meta_key=meta_key)

    @countTime
    def assign_cell_indexes(self,LC):
        for b in self:
            b.assign_cell_index(LC)

    def validate_images(self,box):
        return all([b.validate_image(box) for b in self])

    def get(self,**crits):
        result=[]
        for e in self:
            if e.match(**crits):
                result.append(e)
        return type(self)(result)
    
class PSFAtom(PSFTopoElement):
    def __init__(self,atomline):
        tokens=[x.strip() for x in atomline.split()]
        assert len(tokens)==9,f'cannot parse psf atomline: {atomline}'
        idx=int(tokens[0])
        super().__init__([idx])
        r,i=split_ri(tokens[2])
        self.serial=idx
        self.chainID=tokens[1]
        self.resseqnum=r
        self.insertion=i
        self.resname=tokens[3]
        self.name=tokens[4]
        self.type=tokens[5]
        self.charge=float(tokens[6])
        self.atomicwt=float(tokens[7])
        if len(segtype_of_resname)==0:
            c=Config()
        self.segtype=segtype_of_resname[self.resname]
        self.ligands=[]
    
    def __hash__(self):
        return self.serial
    
    def __eq__(self,other):
        return self.serial==other.serial

    def __str__(self):
        return f'{self.chainID}_{self.resname}_{self.resseqnum}{self.insertion}-{self.name}({self.serial})'

    def isH(self):
        return self.atomicwt<1.1
    
    def add_ligand(self,other):
        if not other in self.ligands:
            self.ligands.append(other)

    def is_pep(self,other):
        return (self.type=='C' and other.type=='NH1') or (self.type=='NH1' and other.type=='C')

class PSFAtomList(PSFTopoElementList):

    def __init__(self,atoms):
        super().__init__(atoms)

    def graph(self):
        g=nx.Graph()
        for a in self:
            for l in a.ligands:
                g.add_edge(a,l)
        logger.debug(f'atomlist graph has {len(g)} nodes')
        return g

class PSFBond(PSFTopoElement):

    def __eq__(self,other):
        return [self.serial1,self.serial2]==[other.serial1,other.serial2] or [self.serial1,self.serial2]==[other.serial2,other.serial1]
    
class PSFBondList(PSFTopoElementList):
    @singledispatchmethod
    def __init__(self,input_data,**kwargs):
        super().__init__(input_data)

    @__init__.register(PSFTopoElementList)
    def _from_tel(self,input_obj,**kwargs):
        super().__init__(input_obj)

    @__init__.register(LineList)
    def _from_list(self,lines,include_serials=[]):
        B=[]
        for line in lines:
            tokens=[int(x) for x in line.split()]
            assert len(tokens)%2==0,f'Poorly formatted BOND line in psffile?'
            if not include_serials:
                for l,r in zip(tokens[:-1:2],tokens[1::2]):
                    B.append(PSFBond([l,r]))
            else:
                for l,r in zip(tokens[:-1:2],tokens[1::2]):
                    if include_serials[l-1] and include_serials[r-1]:
                        B.append(PSFBond([l,r]))
        super().__init__(B)

    def to_graph(self):
        g=nx.Graph()
        for b in self:
            g.add_edge(b.serial1,b.serial2)
        return g

class PSFAngle(PSFTopoElement):

    def __eq__(self,other):
        return [self.serial1,self.serial2,self.serial3]==[other.serial1,other.serial2,other.serial3] or [self.serial1,self.serial2,self.serial3]==[other.serial3,other.serial2,other.serial1] 

class PSFAngleList(PSFTopoElementList):
    @singledispatchmethod
    def __init__(self,input_data,**kwargs):
        super().__init__(input_data)

    @__init__.register(PSFTopoElementList)
    def _from_tel(self,input_obj,**kwargs):
        super().__init__(input_obj)

    @__init__.register(LineList)
    def _from_list(self,lines,include_serials=[]):
        B=[]
        for line in lines:
            tokens=[int(x) for x in line.split()]
            assert len(tokens)%3==0,f'Poorly formatted ANGLE line in psffile?'
            if not include_serials:
                for l,m,r in zip(tokens[:-2:3],tokens[1:-1:3],tokens[2::3]):
                    B.append(PSFAngle([l,m,r]))
            else:
                for l,m,r in zip(tokens[:-2:3],tokens[1:-1:3],tokens[2::3]):
                    if include_serials[l-1] and include_serials[m-1] and include_serials[r-1]:
                        B.append(PSFBond([l,m,r]))
        super().__init__(B)

class PSFDihedral(PSFTopoElement):
    def __eq__(self,other):
        return [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial1,other.serial2,other.serial3,other.serial4] or [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial4,other.serial3,other.serial2,other.serial1] 

class PSFDihedralList(PSFTopoElementList):
    @singledispatchmethod
    def __init__(self,input_data,**kwargs):
        super().__init__(input_data)

    @__init__.register(PSFTopoElementList)
    def _from_tel(self,input_obj,**kwargs):
        super().__init__(input_obj)

    @__init__.register(LineList)
    def _from_list(self,lines,include_serials=[]):
        B=[]
        for line in lines:
            tokens=[int(x) for x in line.split()]
            assert len(tokens)%4==0,f'Poorly formatted PHI line in psffile?'
            if not include_serials:
                for li,i,j,rj in zip(tokens[:-3:4],tokens[1:-2:4],tokens[2:-1:4],tokens[3::4]):
                    B.append(PSFDihedral([li,i,j,rj]))
            else:
                for li,i,j,rj in zip(tokens[:-3:4],tokens[1:-2:4],tokens[2:-1:4],tokens[3::4]):
                    if include_serials[li-1] and include_serials[i-1] and include_serials[j-1] and include_serials[rj-1]:
                        B.append(PSFDihedral([li,i,j,rj]))
        super().__init__(B)
    
class PSFContents:
    def __init__(self,filename,topology_segtypes=[],parse_topology=[]):
        logger.debug(f'Reading {filename}...')
        with open(filename,'r') as f:
            psflines=f.read().split('\n')
        logger.debug(f'{len(psflines)} lines...')
        self.token_idx={}
        self.token_count={}
        self.patches={}
        for i,l in enumerate(psflines):
            toktst=[x.strip() for x in l.split()]
            if len(toktst)>=2 and toktst[1][0]=='!':
                token_name=toktst[1][2:]
                if token_name[-1]==':':
                    token_name=token_name[:-1]
                self.token_idx[token_name]=i
                self.token_count[token_name]=int(toktst[0])
            elif len(toktst)>1 and toktst[0]=="REMARKS":
                remark_type=toktst[1]
                if remark_type=='patch':
                    patch_type=toktst[2]
                    if patch_type not in ['NTER','CTER']:
                        if not patch_type in self.patches:
                            self.patches[patch_type]=[]
                        self.patches[patch_type].append(toktst[3:])
        self.token_lines={}
        for (k,l0),l1 in zip(self.token_idx.items(),list(self.token_idx.values())[1:]+[len(psflines)]):
            fl=l0+1
            ll=l1-1
            self.token_lines[k]=psflines[fl:ll]
        logger.debug(f'{len(self.token_lines)} tokensets:')
        logger.debug(f'{", ".join([x for x in self.token_lines.keys()])}')
        self.atoms=PSFAtomList([PSFAtom(x) for x in self.token_lines['ATOM']])
        self.atomserials=[x.serial for x in self.atoms]
        logger.debug(f'{len(self.atoms)} total atoms...')
        self.ssbonds=SSBondList([SSBond(L) for L in self.patches.get('DISU',[])])
        self.links=LinkList([])
        for patchtype,patchlist in self.patches.items():
            if patchtype in Link.patch_atomnames:
                for patch in patchlist:
                    self.links.append(Link([patchtype]+patch))
        if parse_topology:
            include_serials=[]
            if topology_segtypes:
                self.included_atoms=PSFAtomList([])
                for segtype in topology_segtypes:
                    sublist=self.atoms.get(segtype=segtype)
                    logger.debug(f'{len(sublist)} atoms of segtype {segtype}...')
                    self.included_atoms.extend(sublist)
                include_serials=[x.segtype in topology_segtypes for x in self.atoms]
                logger.debug(f'Including {include_serials.count(True)}/{len(include_serials)} topologically active atoms ({len(self.included_atoms)}) from segtypes {topology_segtypes}')
            if 'bonds' in parse_topology:
                self.bonds=PSFBondList(LineList(self.token_lines['BOND']),include_serials=include_serials)
                # logger.debug(f'Creating graph from {len(self.bonds)} bonds...')
                self.G=self.bonds.to_graph()
                logger.debug(f'extending atom instances with ligands...')
                self.add_ligands()
                logger.debug(f'done')
                # logger.debug(f'Parsed {len(self.bonds)} bonds.')
            if 'angles' in parse_topology:
                self.angles=PSFAngleList(LineList(self.token_lines['THETA']),include_serials=include_serials)
            if 'dihedrals' in parse_topology:
                self.dihedrals=PSFDihedralList(LineList(self.token_lines['PHI']),include_serials=include_serials)
            if 'impropers' in parse_topology:
                self.dihedrals=PSFDihedralList(LineList(self.token_lines['IMPHI']),include_serials=include_serials)
    
    def add_ligands(self):
        assert hasattr(self,'bonds')
        for i,a in enumerate(self.atoms):
            a.ligands=[]
            assert a.serial==i+1
        for b in self.bonds:
            i,j=b.serial1,b.serial2
            aix,ajx=i-1,j-1#self.atomserials.index(i),self.atomserials.index(j)
            ai,aj=self.atoms[aix],self.atoms[ajx]
            ai.add_ligand(aj)
            aj.add_ligand(ai)

    def get_charge(self):
        return np.sum([x.charge for x in self.atoms])