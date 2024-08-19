# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import networkx as nx
from collections import UserDict
from rdkit import Chem
import yaml
import os
import glob
from .config import ResourceManager

import logging
logger=logging.getLogger(__name__)

class CharmmMassRecord:
    def __init__(self,record):
        tok=record
        self.com=''
        tokcom=tok.split('!')
        if len(tokcom)>1:
            tok=tokcom[0]
            self.com='!'.join(tokcom[1:])
        tokens=[t.strip() for t in tok.split()]
        self.atom_type=tokens[2]
        self.atom_mass=float(tokens[3])
        try:
            self.atom_element=tokens[4]
        except:
            self.atom_element=self.atom_type[0]
        # logger.debug(f'{self.atom_type} {self.atom_mass} {self.atom_element}')

class CharmmMasses(UserDict):
    def __init__(self,mList):
        self.data={}
        for m in mList:
            self.data[m.atom_type]=m

class CharmmTopAtom:
    def __init__(self,atomstring,masses=[]):
        tokens=atomstring.split()
        self.name=tokens[1]
        self.type=tokens[2]
        self.charge=float(tokens[3])
        self.mass=0.0
        self.element='?'
        if masses:
            atom_massrecord=masses.get(self.type,None)
            if atom_massrecord:
                self.mass=atom_massrecord.atom_mass
                self.element=atom_massrecord.atom_element
        
class CharmmBond:
    def __init__(self,name1,name2,degree=1):
        self.name1=name1
        self.name2=name2
        self.degree=degree

class CharmmAtoms(UserDict):
    def __init__(self,a_list):
        self.data={a.name:a for a in a_list}

class CharmmTopResi:
    def __init__(self,blockstring,masses=[]):
        lines=blockstring.split('\n')
        titlecard=lines[0]
        tctokens=titlecard.split()
        self.resname=tctokens[1]
        # logger.debug(f'{len(blockstring)} bytes in {len(lines)} lines in resi {self.resname}')
        if tctokens[2].endswith('!'):
            self.charge=float(tctokens[2][:-1])
            tctokens.insert(3,'!')
        else:
            self.charge=float(tctokens[2])
        self.synonym=''
        if '!' in tctokens:
            idx=tctokens.index('!')
            self.synonym=' '.join(tctokens[idx+1:])
        didx=1
        while lines[didx].startswith('!'): didx+=1
        datacards=lines[didx:]
        # logger.debug(f'{len(datacards)} datacards in resi {self.resname}')
        self.groupcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('GROUP')]
        self.bondcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('BOND')]
        self.doublecard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('DOUBLE')]
        # logger.debug(f'groupcards at {self.groupcard_idx}')
        # logger.debug(f'bondcards at {self.bondcard_idx}')
        # logger.debug(f'doublecards at {self.doublecard_idx}')
        self.atoms=[]
        self.atoms_in_group={}
        self.comment_strings=[]
        for gn,gi in enumerate(self.groupcard_idx):
            if '!' in datacards[gi]:
                bang_idx=datacards[gi].index('!')
                comment=datacards[gi][bang_idx+1:]
                if len(comment.strip())>0:
                    self.comment_strings.append(comment)
            ai=1
            while datacards[gi+ai].startswith('ATOM'):
                if '!' in datacards[gi+ai]:
                    bang_idx=datacards[gi+ai].index('!')
                    atom_string=datacards[gi+ai][:bang_idx]
                    comment=datacards[gi+ai][bang_idx+1:]
                else:
                    atom_string=datacards[gi+ai]
                    comment=''
                if not gn in self.atoms_in_group:
                    self.atoms_in_group[gn]=[]
                at=CharmmTopAtom(atom_string,masses=masses)
                self.atoms_in_group[gn].append(at)
                self.atoms.append(at)
                if len(comment.strip())>0:
                    self.comment_strings.append(comment)
                ai+=1
        self.atomdict=CharmmAtoms(self.atoms)
        self.bonds=[]
        for bn,bi in enumerate(self.bondcard_idx):
            bondcard=datacards[bi]
            tokens=bondcard.split()[1:]
            for n1,n2 in zip(tokens[:-1:2],tokens[1::2]):
                b=CharmmBond(n1,n2)
                self.bonds.append(b)
        for bn,bi in enumerate(self.doublecard_idx):
            doublecard=datacards[bi]
            tokens=doublecard.split()[1:]
            for n1,n2 in zip(tokens[:-1:2],tokens[1::2]):
                b=CharmmBond(n1,n2,degree=2)
                self.bonds.append(b)

    def num_atoms(self):
        return len(self.atoms)

    def formula(self):
        fdict={}
        sortorder='CHNOP'
        for a in self.atoms:
            if a.element!='?':
                if not a.element in fdict:
                    fdict[a.element]=0
                fdict[a.element]+=1
        if fdict:
            retstr=''
            for e in sortorder:
                if e in fdict:
                    n='' if fdict[e]==1 else str(fdict[e])
                    retstr+=f'{e}{n}'
            for k,v in fdict.items():
                if not k in sortorder:
                    n='' if v==1 else str(v)
                    retstr+=f'{k}{n}'
        return retstr

    def mass(self):
        sum=0.0
        for a in self.atoms:
            sum+=a.mass
        return sum

    def to_graph(self,includeH=True):
        g=nx.Graph()
        for b in self.bonds:
            if b.name1.startswith('H') or b.name2.startswith('H'):
                if includeH:
                    g.add_edge(b.name1,b.name2)
            else:
                g.add_edge(b.name1,b.name2)
        return g
    
def getMasses(topfile):
    R=[]
    with open(topfile,'r') as f:
        lines=f.read().split('\n')
    masses=[]
    for i in range(len(lines)):
        if lines[i].startswith('MASS'):
            masses.append(CharmmMassRecord(lines[i]))
    return masses

def getResis(topfile,masses=[]):
    R=[]
    with open(topfile,'r') as f:
        lines=f.read().split('\n')
    residx=[]
    for i in range(len(lines)):
        if lines[i].startswith('RESI'):
            residx.append(i)
    if not residx:
        return R
    # logger.debug(f'{topfile}: residx {residx}')
    bufs=[]
    for ll,lr in zip(residx[:-1],residx[1:]):
        bufs.append('\n'.join(lines[ll:lr]))
    bufs.append('\n'.join(lines[residx[-1]:]))
    # logger.debug(f'topfile: {len(bufs)} resis')
    # for i in range(len(bufs)):
    #     logger.debug(f'{len(bufs[i])} bytes in resi {i} {bufs[i][:35]}')
    for block in bufs:
        resi=CharmmTopResi(block,masses=masses)
        resi.topfile=topfile
        R.append(resi)
    return R

def makeBondGraph(mol):
    g2=nx.Graph()
    atoms=[]
    for a in mol.GetAtoms():
        atoms.append(a.GetIdx())
    # logger.debug(f'{len(atoms)} atoms')
    for bond in mol.GetBonds():
        aid1=atoms[bond.GetBeginAtomIdx()]
        aid2=atoms[bond.GetEndAtomIdx()]
        g2.add_edge(aid1,aid2)
    # logger.debug(f'returning {g2}')
    return g2,atoms


# if __name__=='__main__':

#     r=ResourceManager()
#     toppardir=r['charmmff']['toppar']
#     lipidstreamdir=os.path.join(toppardir,'stream/lipid')
#     tops=os.path.join(toppardir,'top_all36_lipid.rtf')+glob.glob(os.path.join(lipidstreamdir,r'*.rtf'))

#     resis=[]
#     for t in tops:
#         resis.extend(getResis(t))

#     for r in resis:
#         print(r.resname,r.synonym,r.num_atoms())

    # with(open('lmsd_charmm.yaml')) as f:
    #     parms=yaml.safe_load(f)
    # # print(parms)
    # toppar=parms['charmmff_dir']
    # sdf=parms['lmsd_database']
    # lipid_db=Chem.SDMolSupplier(sdf)
    # lmids=[a.GetProp('LM_ID') for a in lipid_db]

    # resi_db={}
    # for mapping in parms['mappings']:
    #     idx=lmids.index(mapping['lmsd_id'])
    #     print(idx)
    #     m_heavy=lipid_db[idx]if block.startswith('RESI'):
            
    #     m=Chem.AddHs(m_heavy)
    #     for prop in m.GetPropNames():
    #         print(f'{prop:>12s}: {m.GetProp(prop)}')
    #     if not mapping['charmm_topfile'] in resi_db:
    #         resi_db[mapping['charmm_topfile']]=getResis(os.path.join(toppar,mapping['charmm_topfile']))
    #     resnames=[x.resname for x in resi_db[mapping['charmm_topfile']]]
    #     idx=resnames.index(mapping['charmm_resname'])
    #     top_mol=resi_db[mapping['charmm_topfile']][idx]
    #     charmm_G=top_mol.to_graph()
    #     lmsd_G,atoms=makeBondGraph(m)
    #     rdict=nx.vf2pp_isomorphism(lmsd_G,charmm_G)
    #     if rdict:
    #         new_names=[rdict[x] for x in atoms]
    #         for a,n in zip(m.GetAtoms(),new_names):
    #             info=Chem.rdchem.AtomMonomerInfo()
    #             info.SetName(n)
    #             a.SetMonomerInfo(info)
    #             resinfo=Chem.rdchem.AtomPDBResidueInfo()
    #             resinfo.SetResidueName(mapping['charmm_resname'])
    #             resinfo.SetResidueNumber(1)
    #             a.SetPDBResidueInfo(resinfo)
    #         w=Chem.PDBWriter(f'{mapping["charmm_resname"]}-renamed.pdb')
    #         w.write(m)
    #         w.flush()
    #         w.close()
