# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""PSF file processing
"""
from .basemod import *
from .mods import SSBond,SSBondList,Link,LinkList
from functools import singledispatchmethod
from .stringthings import split_ri
from .config import segtype_of_resname, Config
import networkx as nx 
import numpy as np
import logging

logger=logging.getLogger(__name__)

class PSFAtom(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['serial','chainID','resseqnum','insertion','resname','name','type','charge','atomicwt']
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(str)
    def _from_psfatomline(self,atomline):
        tokens=[x.strip() for x in atomline.split()]
        assert len(tokens)==9,f'cannot parse psf atomline: {atomline}'
        r,i=split_ri(tokens[2])
        input_dict={
            'serial':int(tokens[0]),
            'chainID':tokens[1],
            'resseqnum':r,
            'insertion':i,
            'resname':tokens[3],
            'name':tokens[4],
            'type':tokens[5],
            'charge':float(tokens[6]),
            'atomicwt':float(tokens[7])
        }
        self.ligands=[]
        super().__init__(input_dict)
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

    def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ'],meta_key=[],box=None):
        self.r=A[A[idx_key]==self.serial][pos_key].to_numpy()[0]
        # logger.debug(f'atom r {self.r}')
        self.meta={}
        for key in meta_key:
            self.meta[key]=A[A[idx_key]==self.serial][key].values[0]
            # logger.debug(f'{key} {self.meta[key]}')

class PSFAtomList(AncestorAwareModList):
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(list)
    def _from_psflines(self,lines):
        A=[]
        for line in lines:
            A.append(PSFAtom(line))
        super().__init__(A)

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
    
    def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ'],meta_key=[],box=None):
        for at in self:
            at.injest_coordinates(A,idx_key=idx_key,pos_key=pos_key,meta_key=meta_key,box=box)

class PSFBond(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['serial1','serial2']
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(list)
    def _from_psflines(self,pair):
        input_dict={
            'serial1':pair[0],
            'serial2':pair[1]
        }
        super().__init__(input_dict)
        self.is_attachment_site=False

    def __eq__(self,other):
        return [self.serial1,self.serial2]==[other.serial1,other.serial2] or [self.serial1,self.serial2]==[other.serial2,other.serial1] 
    
    def is_pep(self):
        if self.atom1.segtype=='protein' and self.atom2.segtype=='protein':
            n12=[self.atom1.name,self.atom2.name]
            if n12==['C','N'] or n12==['N','C']:
                return True
        return False
    
    def validate_image(self,box):
        minboxlen=min([box[i][i] for i in [0,1,2]])
        ri=self.atom1.r
        rj=self.atom2.r
        self.COM=0.5*(ri+rj)
        in_same_img=[np.linalg.norm(ri-self.COM)<(minboxlen/2),np.linalg.norm(rj-self.COM)<(minboxlen/2)]
        return all(in_same_img)
    

class PSFBondList(AncestorAwareModList):

    @singledispatchmethod
    def __init__(self,input_obj,**kwargs):
        super().__init__(input_obj)

    @__init__.register(list)
    def _from_psflines(self,lines,include_only=PSFAtomList([])):
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
                        b.atom1=include_only[li]
                        b.atom2=include_only[ri]
                        b.atom1.add_ligand(b.atom2)
                        b.atom2.add_ligand(b.atom1)
                        B.append(b)
                        if len(B)%1000==0:
                            logger.debug(f'{len(B)}...')
                    except:
                        pass
        super().__init__(B)

    def to_graph(self):
        g=nx.Graph()
        for b in self:
            g.add_edge(b.serial1,b.serial2)
        return g

    def validate_images(self,box):
        res=[]
        for b in self:
            res.append(b.validate_image(box))
        return all(res)

class PSFAngle(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['serial1','serial2','serial3']
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(list)
    def _from_psflines(self,triplet):
        input_dict={
            'serial1':triplet[0],
            'serial2':triplet[1],
            'serial3':triplet[2]
        }
        super().__init__(input_dict)

    def __eq__(self,other):
        return [self.serial1,self.serial2,self.serial3]==[other.serial1,other.serial2,other.serial3] or [self.serial1,self.serial2,self.serial3]==[other.serial3,other.serial2,other.serial1] 
        
    def validate_image(self,box):
        minboxlen=min([box[i][i] for i in [0,1,2]])
        ri=self.atom1.r
        rj=self.atom2.r
        rk=self.atom3.r
        self.COM=0.5*(ri+rj+rk)
        in_same_img=[np.linalg.norm(ri-self.COM)<(minboxlen/2),np.linalg.norm(rj-self.COM)<(minboxlen/2),np.linalg.norm(rk-self.COM)<(minboxlen/2)]
        return all(in_same_img)

class PSFAngleList(AncestorAwareModList):

    @singledispatchmethod
    def __init__(self,input_obj,**kwargs):
        super().__init__(input_obj)

    @__init__.register(list)
    def _from_psflines(self,lines,include_only=PSFAtomList([])):
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
                        a.atom1=include_only[ok_serials.index(l)]
                        a.atom2=include_only[ok_serials.index(c)]
                        a.atom3=include_only[ok_serials.index(r)]
                        A.append(a)
                else:
                    A.append(PSFAngle([l,c,r]))

        super().__init__(A)

    def validate_images(self,box):
        res=[]
        for a in self:
            res.append(a.validate_image(box))
        return all(res)

class PSFDihedral(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['serial1','serial2','serial3','serial4']
    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(list)
    def _from_psflines(self,quad):
        input_dict={
            'serial1':quad[0],
            'serial2':quad[1],
            'serial3':quad[2],
            'serial4':quad[3]
        }
        super().__init__(input_dict)

    def __eq__(self,other):
        return [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial1,other.serial2,other.serial3,other.serial4] or [self.serial1,self.serial2,self.serial3,self.serial4]==[other.serial4,other.serial3,other.serial2,other.serial1] 
        
    def validate_image(self,box):
        minboxlen=min([box[i][i] for i in [0,1,2]])
        ri=self.atom1.r
        rj=self.atom2.r
        rk=self.atom3.r
        rl=self.atom4.r
        self.COM=0.5*(ri+rj+rk+rl)
        in_same_img=[np.linalg.norm(ri-self.COM)<(minboxlen/2),np.linalg.norm(rj-self.COM)<(minboxlen/2),np.linalg.norm(rk-self.COM)<(minboxlen/2),np.linalg.norm(rl-self.COM)<(minboxlen/2)]
        return all(in_same_img)

class PSFDihedralList(AncestorAwareModList):

    @singledispatchmethod
    def __init__(self,input_obj,**kwargs):
        super().__init__(input_obj)

    @__init__.register(list)
    def _from_psflines(self,lines,include_only=PSFAtomList([])):
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
                        d.atom1=include_only[ok_serials.index(i)]
                        d.atom2=include_only[ok_serials.index(j)]
                        d.atom3=include_only[ok_serials.index(k)]
                        d.atom4=include_only[ok_serials.index(l)]
                        D.append(d)
                else:
                    D.append(PSFDihedral([i,j,k,l]))

        super().__init__(D)

    def validate_images(self,box):
        res=[]
        for a in self:
            res.append(a.validate_image(box))
        return all(res)
    
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
        self.atoms=PSFAtomList(self.token_lines['ATOM'])
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
            