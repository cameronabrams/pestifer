# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""PSF file processing
"""
from .basemod import *
from .mods import SSBond,SSBondList,Link,LinkList
from .coord import mic_shift
from .stringthings import split_ri
from .config import segtype_of_resname, Config
import networkx as nx 
import numpy as np
import pandas as pd
import logging
from copy import deepcopy
from argparse import Namespace

logger=logging.getLogger(__name__)

class PSFAtom(Namespace):
    def __init__(self,atomline):
        tokens=[x.strip() for x in atomline.split()]
        assert len(tokens)==9,f'cannot parse psf atomline: {atomline}'
        r,i=split_ri(tokens[2])
        super().__init__(
            serial=int(tokens[0]),
            chainID=tokens[1],
            resseqnum=r,
            insertion=i,
            resname=tokens[3],
            name=tokens[4],
            type=tokens[5],
            charge=float(tokens[6]),
            atomicwt=float(tokens[7])
        )
        self.ligands=[]
        if len(segtype_of_resname)==0:
            c=Config()
        self.segtype=segtype_of_resname[self.resname]
        # logger.debug(f'Assigned segtype {self.segtype} to PSFAtom serial {self.serial}')
    
    def __hash__(self):
        return self.serial
    
    def __eq__(self,other):
        return self.serial==other.serial

    def __str__(self):
        return f'{self.chainID}_{self.resname}{self.resseqnum}{self.insertion}-{self.name}({self.serial})'

    def isH(self):
        return self.atomicwt<1.1

    def is_pep(self,other):
        if not other in self.ligands: return False
        if self.segtype=='protein' and other.segtype=='protein':
            n12=[self.name,other.name]
            if n12==['C','N'] or n12==['N','C']:
                return True
        return False

    def add_ligand(self,other):
        if not other in self.ligands:
            self.ligands.append(other)

    def add_attr(self,attrname,attrval):
        self.__dict__[attrname]=attrval

    def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ'],meta_key=[]):
        if idx_key==None:
            self.r=A.loc[self.serial][pos_key].to_numpy()[0]
        else:
            self.r=A[A[idx_key]==self.serial][pos_key].to_numpy()[0]
        # logger.debug(f'atom r {self.r}')
        self.meta={}
        for key in meta_key:
            if idx_key==None:
                self.meta[key]=A.loc[self.serial][key].values[0]
            else:
                self.meta[key]=A[A[idx_key]==self.serial][key].values[0]
            # logger.debug(f'{key} {self.meta[key]}')

    def match(self,**crits):
        return all([self.__dict__[k]==v for k,v in crits.items()])

class PSFAtomList(UserList):

    def __init__(self,atoms):
        super().__init__(atoms)
        self.df=pd.DataFrame({'chainID':[x.chainID for x in self],
                              'resseqnum':[x.resseqnum for x in self],
                              'insertion':[x.insertion for x in self],
                              'resname':[x.resname for x in self],
                              'name':[x.name for x in self],
                              'type':[x.type for x in self],
                              'charge':[x.charge for x in self],
                              'atomicwt':[x.atomicwt for x in self]},
                              index=[x.serial for x in self])
        
    def add_attr(self,attrname,attrval):
        for item in self:
            item.add_attr(attrname,attrval)

    def graph(self):
        g=nx.Graph()
        my_serials=[x.serial for x in self]
        for a in self:
            for l in a.ligands:
                if l.serial in my_serials:
                    g.add_edge(a,l)
        return g
    
    def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ']):
        if idx_key==None:
            self.df[pos_key]=A[pos_key]
            for at,(i,row) in zip(self,self.df.iterrows()):
                at.r=row[pos_key].to_numpy()[0]
        else:
            for at in self:
                at.injest_coordinates(A,idx_key=idx_key,pos_key=pos_key)

    def get(self,**crits):
        result=[]
        for at in self:
            if at.match(**crits):
                result.append(at)
        return type(self)(result)

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

    def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ']):
        if idx_key==None:
            self.P=np.array(A.loc[self.idx_list][pos_key].values)
        else:
            self.P=np.array(A[A[idx_key].isin(self.idx_list)][pos_key].values)
        self.COM=np.mean(self.P,axis=0)

    def validate_image(self,box):
        boxlen=np.array([box[i][i] for i in [0,1,2]])
        in_same_img=np.linalg.norm(self.P-self.COM)<(boxlen/2)
        self.same_image=np.all(in_same_img)
        return self.same_image
    
    def mic_shift(self,reference_point,box):
        for i in range(len(self.P)):
            self.P[i]=mic_shift(self.P[i],reference_point,box)
        return self

    def copy(self):
        return deepcopy(self)
    
class PSFTopoElementList(UserList):
    def __init__(self,data):
        super().__init__(data)

    def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ']):
        for b in self:
            b.injest_coordinates(A,idx_key=idx_key,pos_key=pos_key)

    def validate_images(self,box):
        return all([b.validate_image(box) for b in self])
    
class PSFBond(PSFTopoElement):

    def __eq__(self,other):
        return [self.serial1,self.serial2]==[other.serial1,other.serial2] or [self.serial1,self.serial2]==[other.serial2,other.serial1] 
    
class PSFBondList(PSFTopoElementList):

    def __init__(self,lines,include_only=PSFAtomList([])):
        ok_serials=[]
        if len(include_only)>0:
            ok_serials=[x.serial for x in include_only]
            logger.debug(f'only including from among {len(ok_serials)} atom serials...')
        B=[]
        for line in lines:
            tokens=[int(x) for x in line.split()]
            assert len(tokens)%2==0,f'Poorly formatted BOND line in psffile?'
            for l,r in zip(tokens[:-1:2],tokens[1::2]):
                if not ok_serials: B.append(PSFBond([l,r]))
                else:
                    try:
                        li=ok_serials.index(l)
                        ri=ok_serials.index(r)
                        b=PSFBond([l,r])
                        B.append(b)
                        if len(B)%1000==0:
                            logger.debug(f'{len(B)} psf bonds processed...')
                    except:
                        pass
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

    def __init__(self,lines,include_only=PSFAtomList([])):
        ok_serials=[]
        if len(include_only)>0:
            ok_serials=[x.serial for x in include_only]
            logger.debug(f'only including from among {len(ok_serials)} atom serials...')
        A=[]
        for line in lines:
            tokens=[int(x) for x in line.split()]
            assert len(tokens)%3==0,f'Poorly formatted THETA line in psffile?'
            for l,c,r in zip(tokens[:-2:3],tokens[1:-1:3],tokens[2::3]):
                if len(ok_serials)>0:
                    if all([_ in ok_serials for _ in [l,c,r]]):
                        # logger.debug(f'admitting {l}-{r}...')
                        a=PSFAngle([l,c,r])
                        A.append(a)
                else:
                    A.append(PSFAngle([l,c,r]))
            if len(A)%1000==0:
                logger.debug(f'{len(A)} psf angles processed...')
        super().__init__(A)

class PSFDihedral(PSFTopoElement):
    def __eq__(self,other):
        return [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial1,other.serial2,other.serial3,other.serial4] or [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial4,other.serial3,other.serial2,other.serial1] 

class PSFDihedralList(PSFTopoElementList):
    def __init__(self,lines,include_only=PSFAtomList([])):
        ok_serials=[]
        if len(include_only)>0:
            ok_serials=[x.serial for x in include_only]
            logger.debug(f'only including from among {len(ok_serials)} atom serials...')
        D=[]
        for line in lines:
            tokens=[int(x) for x in line.split()]
            assert len(tokens)%4==0,f'Poorly formatted PHI line in psffile?'
            for i,j,k,l in zip(tokens[:-3:4],tokens[1:-2:4],tokens[2:-1:4],tokens[3::4]):
                if len(ok_serials)>0:
                    if all([_ in ok_serials for _ in [i,j,k,l]]):
                        # logger.debug(f'admitting {l}-{r}...')
                        d=PSFDihedral([i,j,k,l])
                        D.append(d)
                else:
                    D.append(PSFDihedral([i,j,k,l]))
            if len(D)%1000==0:
                logger.debug(f'{len(D)} psf dihedrals processed...')
        super().__init__(D)
    
class PSFContents:
    def __init__(self,filename,include_solvent=False,parse_topology=[]):
        logger.debug(f'Reading {filename}...')
        with open(filename,'r') as f:
            psflines=f.read().split('\n')
        logger.debug(f'{len(psflines)} lines...')
        self.token_idx={}
        self.patches={}
        for i,l in enumerate(psflines):
            toktst=[x.strip() for x in l.split()]
            if len(toktst)>=2 and toktst[1][0]=='!':
                token_name=toktst[1][2:]
                if token_name[-1]==':':
                    token_name=token_name[:-1]
                self.token_idx[token_name]=i
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
        logger.debug(f'{len(self.atoms)} total atoms...')
        self.ssbonds=SSBondList([SSBond(L) for L in self.patches.get('DISU',[])])
        self.links=LinkList([])
        for patchtype,patchlist in self.patches.items():
            if patchtype in Link.patch_atomnames:
                for patch in patchlist:
                    self.links.append(Link([patchtype]+patch))
        if parse_topology:
            logger.debug(f'Parsing all topology information in {filename}...This may take a while.')
            self.protein_atoms=self.atoms.get(segtype='protein')
            self.glycan_atoms=self.atoms.get(segtype='glycan')
            self.lipid_atoms=self.atoms.get(segtype='lipid')
            self.nonsolvent_atoms=PSFAtomList([])
            if self.protein_atoms: self.nonsolvent_atoms.extend(self.protein_atoms)
            if self.glycan_atoms: self.nonsolvent_atoms.extend(self.glycan_atoms)
            if self.lipid_atoms: self.nonsolvent_atoms.extend(self.lipid_atoms)
            logger.debug(f'{len(self.nonsolvent_atoms)} non-solvent atoms')
            if 'bonds' in parse_topology:
                self.bonds=PSFBondList(self.token_lines['BOND'],include_only=(self.nonsolvent_atoms if not include_solvent else []))
                logger.debug(f'Creating graph from {len(self.bonds)} bonds...')
                self.G=self.bonds.to_graph()
                logger.debug(f'Parsed {len(self.bonds)} bonds.')
            if 'angles' in parse_topology:
                self.angles=PSFAngleList(self.token_lines['THETA'],include_only=(self.nonsolvent_atoms if not include_solvent else []))
            if 'dihedrals' in parse_topology:
                self.dihedrals=PSFDihedralList(self.token_lines['PHI'],include_only=(self.nonsolvent_atoms if not include_solvent else []))
            if 'impropers' in parse_topology:
                self.dihedrals=PSFDihedralList(self.token_lines['IMPHI'],include_only=(self.nonsolvent_atoms if not include_solvent else []))
            