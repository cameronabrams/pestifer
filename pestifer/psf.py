# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""PSF file processing
"""
from .basemod import *
from .mods import SSBond,SSBondList,Link,LinkList
from functools import singledispatchmethod
from .stringthings import split_ri
from .config import segtype_of_resname, Config
import networkx as nx 

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

    def graph(self):
        g=nx.Graph()
        my_serials=[x.serial for x in self]
        for a in self:
            for l in a.ligands:
                if l.serial in my_serials:
                    g.add_edge(a,l)
        return g

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
                if len(ok_serials)>0:
                    if l in ok_serials and r in ok_serials:
                        # logger.debug(f'admitting {l}-{r}...')
                        b=PSFBond([l,r])
                        b.atom1=include_only[ok_serials.index(l)]
                        b.atom2=include_only[ok_serials.index(r)]
                        b.atom1.add_ligand(b.atom2)
                        b.atom2.add_ligand(b.atom1)
                        B.append(b)
                else:
                    B.append(PSFBond([l,r]))

        super().__init__(B)

class PSFContents:
    def __init__(self,filename,include_solvent=False,parse_topology=False):
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
            nonsolvent_atoms=self.atoms.get(segtype='protein')+self.atoms.get(segtype='glycan')
            self.bonds=PSFBondList(self.token_lines['BOND'],include_only=(nonsolvent_atoms if not include_solvent else []))
            logger.debug(f'Parsed {len(self.bonds)} bonds...')
            # self.angles=PSFAngleList(self.token_lines['THETA'],include_only=(nonsolvent_atoms if not include_solvent else []))
            # self.dihedrals=PSFDihedralList(self.token_lines['PHI'],include_only=(nonsolvent_atoms if not include_solvent else []))
            