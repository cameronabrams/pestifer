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

    def injest_coordinates(self,A,pos_key=['posX','posY','posZ'],meta_key=[]):
        self.P=np.array(A.loc[self.idx_list][pos_key].values)
        if len(self.P)>1:
            self.COM=np.mean(self.P,axis=0)
            self.radii=self.P-self.COM
            self.radii/=np.linalg.norm(self.radii)
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
        # shifts=[]
        holder=self.copy()
        for i in range(len(self.P)):
            holder.P[i],bl=mic_shift(self.P[i],reference_point,box)
            # shifts.append(np.any(bl))
        holder.COM=np.mean(holder.P,axis=0)
        holder.radii=holder.P-holder.COM
        holder.radii/=np.linalg.norm(holder.radii)
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
        # self.ligands=[]
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

    # def is_pep(self,other):
    #     if not other in self.ligands: return False
    #     if self.segtype=='protein' and other.segtype=='protein':
    #         n12=[self.name,other.name]
    #         if n12==['C','N'] or n12==['N','C']:
    #             return True
    #     return False

    # def add_ligand(self,other):
    #     if not other in self.ligands:
    #         self.ligands.append(other)

    # def add_attr(self,attrname,attrval):
    #     self.__dict__[attrname]=attrval

    # def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ'],meta_key=[]):
    #     if idx_key==None:
    #         self.r=A.loc[self.serial][pos_key].to_numpy()
    #     else:
    #         self.r=A[A[idx_key]==self.serial][pos_key].to_numpy()[0]
    #     # logger.debug(f'atom r {self.r}')
    #     self.meta={}
    #     for key in meta_key:
    #         if idx_key==None:
    #             self.meta[key]=A.loc[self.serial][key].values[0]
    #         else:
    #             self.meta[key]=A[A[idx_key]==self.serial][key].values[0]
    #         # logger.debug(f'{key} {self.meta[key]}')


class PSFAtomList(PSFTopoElementList):

    def __init__(self,atoms):
        super().__init__(atoms)

    # def make_df(self):
    #     self.df=pd.DataFrame({'chainID':[x.chainID for x in self],
    #                           'resseqnum':[x.resseqnum for x in self],
    #                           'insertion':[x.insertion for x in self],
    #                           'resname':[x.resname for x in self],
    #                           'name':[x.name for x in self],
    #                           'type':[x.type for x in self],
    #                           'charge':[x.charge for x in self],
    #                           'atomicwt':[x.atomicwt for x in self]},
    #                           index=[x.serial for x in self])
        # my_logger(self.df,logger.debug,dfoutmode='info',just='<',fill='')
        # logger.debug(f'\n{self.df.head()}')
        
    # def add_attr(self,attrname,attrval):
    #     for item in self:
    #         item.add_attr(attrname,attrval)

    # def graph(self):
    #     g=nx.Graph()
    #     my_serials=[x.serial for x in self]
    #     for a in self:
    #         for l in a.ligands:
    #             if l.serial in my_serials:
    #                 g.add_edge(a,l)
    #     return g
    
    # @countTime
    # def injest_coordinates(self,A,pos_key=['posX','posY','posZ'],meta_key=[]):
    #     self.df[pos_key]=A[pos_key]
    #     for at,(i,row) in zip(self,A.iterrows()):
    #         at.r=row[pos_key].to_numpy()
    #         for k in meta_key:
    #             at.__dict__[k]=row[k]

    

class PSFBond(PSFTopoElement):

    def __eq__(self,other):
        return [self.serial1,self.serial2]==[other.serial1,other.serial2] or [self.serial1,self.serial2]==[other.serial2,other.serial1]
    
    # def nv(self):
    #     vp=self.P[0]-self.COM
    #     vp/=np.linalg.norm(vp)
    #     vm=self.P[1]-self.COM
    #     vm/=np.linalg.norm(vp)
    #     return vp,vm
    
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
                        # if len(B)%1000==0:
                        #     logger.debug(f'{len(B)} psf bonds processed...')
        super().__init__(B)

    def to_graph(self):
        g=nx.Graph()
        for b in self:
            g.add_edge(b.serial1,b.serial2)
        return g

# class PSFBondDataFrame:
#     @singledispatchmethod
#     def __init__(self,input_obj,**kwargs):
#         pass
#     @__init__.register(list)
#     def from_linelist(self,line_list,include_serials=[]):
#         loc=0
#         nbond_est=4*len(line_list)
#         ii=[0]*nbond_est
#         jj=[0]*nbond_est
#         for line in line_list:
#             tokens=[int(x) for x in line.split()]
#             assert len(tokens)%2==0,f'Poorly formatted BOND line in psffile?'
#             if not include_serials:
#                 for l,r in zip(tokens[:-1:2],tokens[1::2]):
#                     ii[loc]=l
#                     jj[loc]=r
#                     loc+=1
#             else:
#                 for l,r in zip(tokens[:-1:2],tokens[1::2]):
#                     if include_serials[l-1] and include_serials[r-1]:
#                         ii[loc]=l
#                         jj[loc]=r
#                         loc+=1
#         # logger.debug(f'processed {loc} bonds into parallel lists')
#         ii=ii[:loc]
#         jj=jj[:loc]
#         self.df=pd.DataFrame({'i':ii,'j':jj},dtype=int)
#         # logger.debug(f'dataframe contains {self.df.shape[0]} rows')
#         # buf=StringIO()
#         # self.df.info(buf=buf)
#         # logger.debug(buf.getvalue())
#     @__init__.register(pd.DataFrame)
#     def from_dataframe(self,df):
#         self.df=df

#     @countTime
#     def injest_coordinates(self,A,pos_key=['posX','posY','posZ'],meta_key=[],iter=1000):
#         b=A[pos_key]
#         m=A[meta_key]
#         self.df[['ix','iy','iz','jx','jy','jz','comx','comy','comz','nvpx','nvpy','nvpz','nvmx','nvmy','nvmz']]=np.zeros(shape=[self.df.shape[0],15],dtype=float)
#         self.df[meta_key]=[['none']*len(meta_key)]*self.df.shape[0]
#         for i in self.df.index:
#             ai,aj=int(self.df.loc[i,'i']),int(self.df.loc[i,'j'])
#             # ri,rj=b.loc[ai].to_numpy(),b.loc[aj].to_numpy()
#             # self.df.loc[i,['ix','iy','iz']]=ri
#             # self.df.loc[i,['jx','jy','jz']]=rj
#             com=np.mean(np.array([b.loc[ai].to_numpy(),b.loc[aj].to_numpy()]),axis=0)
#             self.df.loc[i,['comx','comy','comz']]=com
#             # vp=ri-com
#             # self.df.loc[i,['nvpx','nvpy','nvpz']]=vp/np.linalg.norm(vp)
#             # vm=rj-com
#             # self.df.loc[i,['nvmx','nvmy','nvmz']]=vm/np.linalg.norm(vm)
#             self.df.loc[i,meta_key]=m.loc[ai]
#             if i>1 and i%iter==0: logger.debug(f'{i} bond coords populated')
#         # logger.debug(f'\n{self.df.head()}')
    
#     def to_graph(self):
#         g=nx.Graph()
#         for i,j in self.df[['i','j']].to_numpy():
#             g.add_edge(i,j)
#         return g

#     @countTime
#     def assign_cell_indexes(self,LC):
#         self.df['linkcell_idx']=[0]*self.df.shape[0]
#         for i in self.df.index:
#             p=self.df.loc[i,['comx','comy','comz']].to_numpy()
#             self.df.loc[i,'linkcell_idx']=LC.ldx_of_cellndx(LC.cellndx_of_point(p))

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
    def __init__(self,filename,topology_segtypes=['protein','glycan','lipid'],parse_topology=[]):
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
        logger.debug(f'{len(self.atoms)} total atoms...')
        self.ssbonds=SSBondList([SSBond(L) for L in self.patches.get('DISU',[])])
        self.links=LinkList([])
        for patchtype,patchlist in self.patches.items():
            if patchtype in Link.patch_atomnames:
                for patch in patchlist:
                    self.links.append(Link([patchtype]+patch))
        if parse_topology:
            # logger.debug(f'Parsing all topology information in {filename}...This may take a while.')
            self.included_atoms=PSFAtomList([])
            for segtype in topology_segtypes:
                self.included_atoms.extend(self.atoms.get(segtype=segtype))
            # self.protein_atoms=self.atoms.get(segtype='protein')
            # self.glycan_atoms=self.atoms.get(segtype='glycan')
            # self.lipid_atoms=self.atoms.get(segtype='lipid')
            # if self.protein_atoms: self.nonsolvent_atoms.extend(self.protein_atoms)
            # if self.glycan_atoms: self.nonsolvent_atoms.extend(self.glycan_atoms)
            # if self.lipid_atoms: self.nonsolvent_atoms.extend(self.lipid_atoms)
            # logger.debug(f'{len(self.included_atoms)} topologically included atoms')
            include_serials=[x.segtype in topology_segtypes for x in self.atoms]
            logger.debug(f'Including {include_serials.count(True)}/{len(include_serials)} topologically active atoms')
            if 'bonds' in parse_topology:
                self.bonds=PSFBondList(LineList(self.token_lines['BOND']),include_serials=include_serials)
                # self.bonddf=PSFBondDataFrame(self.token_lines['BOND'],include_serials=include_serials)
                logger.debug(f'Creating graph from {len(self.bonds)} bonds...')
                self.G=self.bonds.to_graph()
                logger.debug(f'Parsed {len(self.bonds)} bonds.')
            # if 'angles' in parse_topology:
            #     self.angles=PSFAngleList(self.token_lines['THETA'],include_only=(self.nonsolvent_atoms if not include_solvent else []))
            # if 'dihedrals' in parse_topology:
            #     self.dihedrals=PSFDihedralList(self.token_lines['PHI'],include_only=(self.nonsolvent_atoms if not include_solvent else []))
            # if 'impropers' in parse_topology:
            #     self.dihedrals=PSFDihedralList(self.token_lines['IMPHI'],include_only=(self.nonsolvent_atoms if not include_solvent else []))
            