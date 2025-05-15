# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import numpy as np
import os
from .psftopoelement import LineList
from .psfatom import PSFAtom, PSFAtomList
from .psfbond import PSFBondList
from .psfangle import PSFAngleList
from .psfdihedral import PSFDihedralList

from ..objs.ssbond import SSBond, SSBondList
from ..objs.link import Link, LinkList

logger=logging.getLogger(__name__)

def get_toppar_from_psf(filename):
    lines=[]
    with open(filename,'r') as f:
        for line in f:
            if '!NATOM' in line:
                break
            lines.append(line)
    toppars=[]
    for l in lines:
        if l[0:17]==' REMARKS topology':
            tokens=l.split()
            top=tokens[2]
            if 'toppar_' in top:
                fn=os.path.basename(top)
                toppars.append(fn)
    return list(set(toppars))  

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