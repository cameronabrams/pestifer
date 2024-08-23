# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import networkx as nx
from collections import UserDict, UserList
import numpy as np
from .stringthings import my_logger
from .config import ResourceManager
from .controller import Controller
from .stringthings import linesplit
from .scriptwriters import Psfgen,VMD
import glob 
import os
import logging
import shutil
from itertools import combinations
import pandas as pd
from pidibble.pdbparse import PDBParser

logger=logging.getLogger(__name__)

class CharmmMassRecord:
    def __init__(self,record):
        tok=record
        tok,self.com=linesplit(record)
        tokens=[t.strip() for t in tok.split()]
        self.atom_type=tokens[2].upper()
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
        self.name=tokens[1].upper()
        self.type=tokens[2].upper()
        self.charge=float(tokens[3])
        self.mass=0.0
        self.element='?'
        if masses:
            atom_massrecord=masses.get(self.type,None)
            if atom_massrecord is not None:
                self.mass=atom_massrecord.atom_mass
                self.element=atom_massrecord.atom_element
            else:
                logger.debug(f'no mass record for atom type {self.type} (raw {tokens[2]})')
                exit(-1)
        else:
            logger.debug(f'cannot initialize without masses')
            exit(-1)
    
    def __eq__(self,other):
        return self.name==other.name
    
class CharmmTopAtomList(UserList):
    def __init__(self,data):
        self.data=data

    def __next__(self):
        for i in range(len(self)):
            yield self[i]

    def get_atom(self,name):
        L=[x.name for x in self]
        # logger.debug(f'looking for {name} in {L}')
        try:
            idx=L.index(name)
            return self[idx]
        except:
            return None

    def get_serial(self,name):
        a=self.get_atom(name)
        if a is not None:
            if hasattr(a,'serial'):
                return a.serial
        return -1
    
    def get_element(self,name):
        a=self.get_atom(name)
        if a is not None:
            if hasattr(a,'element'):
                if a.element=='?':
                    logger.error(f'{a.name} has no element?')
                # logger.debug(f'returning {a.element} for {name}')
                return a.element
        else:
            logger.debug(f'Could not find atom with name {name}')
        logger.error(f'{name} has no element?')
        return '?'
    
    def append(self,atom):
        if not hasattr(atom,'serial'):
            atom.serial=len(self)+1
        super().append(atom)


class CharmmBond:
    def __init__(self,name1,name2,degree=1):
        self.name1=name1
        self.name2=name2
        self.degree=degree

    def __eq__(self,other):
        c1=self.name1==other.name1 and self.name2==other.name2
        c2=self.name1==other.name2 and self.name2==other.name1
        return c1 or c2
    
class CharmmBondList(UserList):
    def __init__(self,data):
        self.data=data

class CharmmAngle:
    def __init__(self,name1,name2,name3):
        self.name1=name1
        self.name2=name2
        self.name3=name3

    def __eq__(self,other):
        r1=self.name2==other.name2
        c1=self.name1==other.name1 and self.name3==other.name3
        c2=self.name1==other.name3 and self.name3==other.name1
        return r1 and (c1 or c2)

class CharmmAngleList(UserList):
    def __init__(self,data):
        self.data=data

class CharmmDihedral:
    def __init__(self,name1,name2,name3,name4,dihedral_type='proper'):
        self.name1=name1
        self.name2=name2
        self.name3=name3
        self.name4=name4
        self.dihedral_type=dihedral_type

    """ dihedral atom names are ordered in ICs, so we preserve this here """
    def __eq__(self,other):
        l1=[self.name1,self.name2,self.name3,self.name4]
        l2=[other.name1,other.name2,other.name3,other.name4]
        if self.dihedral_type!=other.dihedral_type:
            return False
        return l1==l2

class CharmmDihedralList(UserList):
    def __init__(self,data):
        self.data=data
        
class CharmmTopIC:
    def __init__(self,ICstring):
        ICstring,dummy=linesplit(ICstring)
        toks=[x.strip(' !') for x in ICstring.split()]
        data=np.array(list(map(float,toks[5:])))
        if data[0]==0.0 or data[1]==0.0 or data[3]==0.0 or data[4]==0.0:
            logger.debug(f'{toks[1:5]}: missing ic data: {data}')
            self.empty=True
            return
        self.atoms=toks[1:5]
        if self.atoms[2].startswith('*'):
            self.atoms[2]=self.atoms[2][1:]
            self.dihedral_type='improper'
        else:
            self.dihedral_type='proper'
        if self.dihedral_type=='proper':
            b0=CharmmBond(self.atoms[0],self.atoms[1])
            b0.length=float(toks[5])
            b1=CharmmBond(self.atoms[2],self.atoms[3])
            b1.length=float(toks[9])
            self.bonds=CharmmBondList([b0,b1])
            a1=CharmmAngle(self.atoms[0],self.atoms[1],self.atoms[2])
            a1.degrees=float(toks[6])
            a2=CharmmAngle(self.atoms[1],self.atoms[2],self.atoms[3])
            a2.degrees=float(toks[8])
            self.angles=CharmmAngleList([a1,a2])
            self.dihedral=CharmmDihedral(self.atoms[0],self.atoms[1],self.atoms[2],self.atoms[3])
            self.dihedral.degrees=float(toks[7])
        else:
            b0=CharmmBond(self.atoms[0],self.atoms[2])
            b0.length=float(toks[5])
            b1=CharmmBond(self.atoms[2],self.atoms[3])
            b1.length=float(toks[9])
            self.bonds=CharmmBondList([b0,b1])
            a1=CharmmAngle(self.atoms[0],self.atoms[2],self.atoms[3])
            a1.degrees=float(toks[6])
            a2=CharmmAngle(self.atoms[1],self.atoms[2],self.atoms[3])
            a2.degrees=float(toks[8])
            self.angles=CharmmAngleList([a1,a2])
            self.dihedral=CharmmDihedral(self.atoms[0],self.atoms[1],self.atoms[2],self.atoms[3],dihedral_type='improper')
            self.dihedral.degrees=float(toks[7])
        self.empty=False

    def __str__(self):
        ans=','.join([f'{i+1}{self.atoms[i]}' for i in range(len(self.atoms))])
        ans+=':'
        ans+=','.join([f'{b.length:.4f}' for b in self.bonds])
        ans+=':'
        ans+=','.join([f'{a.degrees:.4f}' for a in self.angles])
        ans+=':'
        ans+=f'{self.dihedral.degrees}-{self.dihedral.dihedral_type}'
        return ans

class CharmmTopResi:
    def __init__(self,blockstring,masses=[]):
        self.blockstring=blockstring
        lines=[x.strip() for x in blockstring.split('\n')]
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
        datacards=[x.upper() for x in lines[didx:]]
        # logger.debug(f'{len(datacards)} datacards in resi {self.resname}')
        self.atomcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('ATOM')]
        self.groupcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('GROU')]
        self.bondcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('BOND')]
        self.doublecard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('DOUB')]
        self.triplecard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('TRIP')]
        self.ICcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('IC')]
        self.atoms=CharmmTopAtomList([])
        self.atoms_in_group={i:CharmmTopAtomList([]) for i in range(len(self.groupcard_idx))}
        self.comment_strings=[]
        if self.atomcard_idx[0]<self.groupcard_idx[0]:
            gn=-1
            logger.debug(f'Orphaned atoms outside of group -- we\'ll just make one')
            ai=self.atomcard_idx[0]
            while datacards[ai].startswith('ATOM') or datacards[ai].startswith('Atom') or datacards[ai].startswith('!') or len(datacards[ai])==0:
                if datacards[ai].startswith('!') or len(datacards[ai])==0:
                    logger.debug(f'non-atom line [{datacards[ai]}]')
                else:
                    atomstring,comment=linesplit(datacards[ai])
                    if not gn in self.atoms_in_group:
                        self.atoms_in_group[gn]=[]
                    at=CharmmTopAtom(atomstring,masses=masses)
                    self.atoms.append(at)
                    self.atoms_in_group[gn].append(at)
                    if len(comment.strip())>0:
                        self.comment_strings.append(comment)
                ai+=1
        if self.groupcard_idx:
            for gn,gi in enumerate(self.groupcard_idx):
                dummy,comment=linesplit(datacards[gi])
                if len(comment.strip())>0:
                    self.comment_strings.append(comment)
                ai=1

                while datacards[gi+ai].startswith('ATOM') or datacards[gi+ai].startswith('Atom') or datacards[gi+ai].startswith('!') or len(datacards[gi+ai])==0:
                    if datacards[gi+ai].startswith('!') or len(datacards[gi+ai])==0:
                        logger.debug(f'non-atom line [{datacards[gi+ai]}]')
                    else:
                        atomstring,comment=linesplit(datacards[gi+ai])
                        if not gn in self.atoms_in_group:
                            self.atoms_in_group[gn]=[]
                        at=CharmmTopAtom(atomstring,masses=masses)
                        self.atoms.append(at)
                        self.atoms_in_group[gn].append(at)
                        if len(comment.strip())>0:
                            self.comment_strings.append(comment)
                    ai+=1
        else:
            for an,ai in enumerate(self.atomcard_idx):
                atomstring,comment=linesplit(datacards[ai])
                at=CharmmTopAtom(atomstring,masses=masses)
                self.atoms.append(at)
                if len(comment.strip())>0:
                    self.comment_strings.append(comment)
                ai+=1
        self.bonds=CharmmBondList([])
        for bn,bi in enumerate(self.bondcard_idx):
            bondcard=datacards[bi].upper()
            tokens=bondcard.split()[1:]
            for n1,n2 in zip(tokens[:-1:2],tokens[1::2]):
                b=CharmmBond(n1,n2)
                self.bonds.append(b)
        for bn,bi in enumerate(self.doublecard_idx):
            doublecard=datacards[bi].upper()
            tokens=doublecard.split()[1:]
            for n1,n2 in zip(tokens[:-1:2],tokens[1::2]):
                b=CharmmBond(n1,n2,degree=2)
                self.bonds.append(b)
        for tn,ti in enumerate(self.triplecard_idx):
            triplecard=datacards[ti].upper()
            tokens=triplecard.split()[1:]
            for n1,n2 in zip(tokens[:-1:2],tokens[1::2]):
                b=CharmmBond(n1,n2,degree=3)
                self.bonds.append(b)
        self.IC=[]
        for icn,ici in enumerate(self.ICcard_idx):
            ICcard=datacards[ici]
            IC=CharmmTopIC(ICcard)
            self.IC.append(IC)

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
            node1=b.name1
            node2=b.name2
            element1=self.atoms.get_element(node1)
            element2=self.atoms.get_element(node2)
            assert element1!='?',f'no element for atom {node1} in {self.resname}'
            assert element2!='?',f'no element for atom {node2} in {self.resname}'
            if element1=='H' or element2=='H':
                if includeH:
                    g.add_edge(node1,node2)
                    g.nodes[node1]['element']=element1
                    g.nodes[node2]['element']=element2
            else:
                g.add_edge(node1,node2)
                g.nodes[node1]['element']=element1
                g.nodes[node2]['element']=element2
        return g
    
    def to_file(self,f):
        f.write(self.blockstring)

    def to_psfgen(self,W):
        # find the first heavy atom that is atom-i in an IC that is not an improper
        refatom=None
        logger.debug(f'trying to find a reference atom from {len(self.atoms)} atoms in this resi')
        i=0
        while refatom is None and i<len(self.atoms):
            atom=self.atoms[i]
            if atom.element!='H':
                for ic in self.IC:
                    if not ic.empty:
                        if ic.atoms[1]==atom.name and ic.dihedral_type != 'improper':
                            refatom=atom
                            refic=ic
            i+=1
        if refatom is None:
            logger.debug(f'Unable to identify a reference atom')
            return -1
        logger.debug(f'refatom is {refatom.name}')
        W.addline(r'segment A {')
        W.addline(r'    first none')
        W.addline(r'    last none')
        W.addline(f'    resid 1 {self.resname}')
        W.addline(r'}')
        # put reference atom at origin
        W.addline(f'psfset coord A 1 {refatom.name} '+r'{0 0 0}')
        b1a1,b1a2,b1l=refic.bonds[0].name1,refic.bonds[0].name2,refic.bonds[0].length
        partner=b1a1
        assert b1a1!=refatom.name
        assert b1a2==refatom.name
        # partner is the first atom in the IC, put along x-axis
        W.addline(f'psfset coord A 1 {partner} '+r'{'+f'{b1l:.5f}'+r' 0 0}')
        a1a1,a1a2,a1a3,a1d=refic.angles[0].name1,refic.angles[0].name2,refic.angles[0].name3,refic.angles[0].degrees
        assert a1a1==partner
        assert a1a2==refatom.name,f'{refatom.name} {refic.atoms}: {refic.angles[0].name1} {refic.angles[0].name2} {refic.angles[0].name3}'
        assert a1a3!=partner
        W.addline(f'set a [expr {a1d}*acos(-1.0)/180.0]')
        W.addline(r'set x [expr 1.5283*cos($a)]')
        W.addline(r'set y [expr 1.5283*sin($a)]')
        W.addline(f'psfset coord A 1 {a1a3} [list $x $y 0.0]')
        return 0

def getMasses(topfile):
    R=[]
    with open(topfile,'r') as f:
        lines=f.read().upper().split('\n')
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
        logger.debug(f'No RESI\'s found in {topfile}')
        return R
    # logger.debug(f'{topfile}: residx {residx}')
    bufs=[]
    for ll,lr in zip(residx[:-1],residx[1:]):
        bufs.append('\n'.join(lines[ll:lr]))
    endidx=len(lines)
    for i,l in enumerate(lines[residx[-1]:]):
        ll=l.upper()
        if ll.startswith('END'):
            oldidx=endidx
            endidx=residx[-1]+i
            # print(topfile,ll,i,residx[-1],endidx,oldidx)
            break
    bufs.append('\n'.join(lines[residx[-1]:endidx]))
    # logger.debug(f'topfile: {len(bufs)} resis')
    # for i in range(len(bufs)):
    #     logger.debug(f'{len(bufs[i])} bytes in resi {i} {bufs[i][:35]}')
    for block in bufs:
        resi=CharmmTopResi(block,masses=masses)
        resi.topfile=topfile
        logger.debug(f'resi {resi.resname} parsed')
        R.append(resi)
    return R

def makeBondGraph_rdkit(mol):
    g2=nx.Graph()
    atoms=[]
    elements=[]
    for a in mol.GetAtoms():
        atoms.append(a.GetIdx())
        elements.append(a.GetSymbol())
    # logger.debug(f'{len(atoms)} atoms')
    for bond in mol.GetBonds():
        aid1=atoms[bond.GetBeginAtomIdx()]
        aid2=atoms[bond.GetEndAtomIdx()]
        g2.add_edge(aid1,aid2)
        g2.nodes[aid1]['element']=elements[bond.GetBeginAtomIdx()]
        g2.nodes[aid2]['element']=elements[bond.GetEndAtomIdx()]
    # logger.debug(f'returning {g2}')
    return g2,atoms


class CharmmResiDatabase(UserDict):
    def __init__(self,streams=[]):
        self.resources=ResourceManager()
        self.toppardir=self.resources['charmmff']['toppar']
        # get all atom types from standard charmmff topology files
        rtftops=glob.glob(os.path.join(self.toppardir,'top_*.rtf'))
        M=[]
        for f in rtftops:
            M.extend(getMasses(f))
        self.all_charmm_topology_files=rtftops
        stream_topology_files={}
        for stream in streams:
            stream_rtfs=glob.glob(os.path.join(self.toppardir,f'*{stream}*.rtf'))
            streamdir=os.path.join(self.toppardir,f'stream/{stream}')
            stream_strs=glob.glob(os.path.join(streamdir,'*.str'))
            stream_strs=[x for x in stream_strs if not ('list' in x or 'model' in x or 'lipid_prot' in x)]
            logger.debug(f'stream strs {stream_strs}')
            stream_topology_files[stream]=stream_rtfs+stream_strs
            self.all_charmm_topology_files.extend(stream_topology_files[stream])
            for f in stream_strs: # already got masses from all rtf files
                M.extend(getMasses(f))
        
        self.M=CharmmMasses(M)
        data={}
        self.charmm_resnames=[]
        nfiles=0
        for stream in streams:
            data[stream]={}
            for t in stream_topology_files[stream]:
                sublist=getResis(t,self.M)
                nfiles+=1
                for resi in sublist:
                    resi.charmmtopfile=t
                    data[stream][resi.resname]=resi
                    self.charmm_resnames.append(resi.resname)

        super().__init__(data)
        logger.debug(f'{nfiles} CHARMM topology files scanned')
        self.charmm_resnames.sort()
        logger.debug(f'{len(self.charmm_resnames)} RESI\'s parsed')
        logger.debug(f'{len(self.M)} atomic mass records parsed')
    
    def get_charmm_topfile(self,charmm_resid):
        if charmm_resid not in self.charmm_resnames:
            return None
        for stream,data in self.items():
            elem=data.get(charmm_resid,None)
            if elem is not None:
                return elem.charmmtopfile
        return None
    
    def get_topo(self,charmm_resid):
        if charmm_resid not in self.charmm_resnames:
            logger.debug(f'resi {charmm_resid} is nowhere')
            return None
        for stream,data in self.items():
            elem=data.get(charmm_resid,None)
            if elem is not None:
                return elem
            else:
                logger.debug(f'resi {charmm_resid} is not in stream {stream}')
        return None

def head_tail_atoms(topo):
    G=topo.to_graph(includeH=False)
    cc=list(nx.chordless_cycles(G))
    logger.debug(f'cc {cc}')
    shortest_paths={}
    heads=[]
    tails=[]
    ring_sizes=[len(list(x)) for x in list(cc)]
    if len(cc)==4 and ring_sizes.count(6)==3 and ring_sizes.count(5)==1: # likely a sterol

        # identify head as carbon-3 of ring A (to which OH is bound)
        # identify tail as carbon-16 on ring D
        ringD=[c for c in cc if len(c)==5][0]
        logger.debug(f'ringD {ringD}')
        edge_partner_idx=[]
        for i in range(len(cc)):
            ci=cc[i]
            logger.debug(f'ring {i} {ci}')
            for j in range(i+1,len(cc)):
                cj=cc[j]
                shared_nodes_count=0
                for ai in ci:
                    for aj in cj:
                        if ai==aj:
                            shared_nodes_count+=1
                if shared_nodes_count==2:
                    edge_partner_idx.append([i,j])
                            # cc[i].remove(ai)
                            # cc[j].remove(ai)
                            # G.remove_node(ai)
        logger.debug(f'edge_partner_index {edge_partner_idx}')
        for i in range(len(cc)):
            nsharededges=0
            for ep in edge_partner_idx:
                if i in ep:
                    nsharededges+=1
            logger.debug(f'ring {i} nsharededges {nsharededges} len {len(cc[i])}')
            if nsharededges==1 and len(cc[i])==6:
                ringA=cc[i]
            elif nsharededges==1 and len(cc[i])==5:
                assert ringD==cc[i]
        logger.debug(f'ringA {ringA}')
        for a in ringA:
            for n in nx.neighbors(G,a):
                if G.nodes[n]['element']=='O':
                    heads=[n]
        paths=[]
        for a in G.__iter__():
            paths.append(len(nx.shortest_path(G,a,heads[0])))
        paths=np.array(paths)
        l_idx=np.argsort(paths)[-1]
        tails=[list(G.__iter__())[l_idx]]
    else:
        ends=[]
        for n in G.__iter__():
            nn=len(list(G.neighbors(n)))
            if nn==1:
                ends.append(n)
        ecount={}
        for e in ends:
            nn=list(G.neighbors(e))[0]
            if not nn in ecount:
                ecount[nn]=[]
            ecount[nn].append(e)
        logger.debug(f'ecount {ecount}')
        for k,v in ecount.items():
            if len(v)>1:
                for vv in v:
                    ends.remove(vv)
                ends.append(k)
        
        logger.debug(f'ends {ends}')
        chain_ends=[]
        DG=G.to_directed()
        for n in ends:
            if G.nodes[n]['element']=='C':
                traversed_edges=list(nx.dfs_edges(DG,source=n,depth_limit=4))
                logger.debug(f'node {n} traversed_edges {traversed_edges}')
                if len(traversed_edges)>=4:
                    carbon_chain=True
                    for e in traversed_edges:
                        a,b=e
                        if G.nodes[a]['element']!='C' or G.nodes[b]['element']!='C':
                            carbon_chain=False
                    if carbon_chain: chain_ends.append(n)
        logger.debug(f'chain_ends: {chain_ends}')
        for e in chain_ends:
            ends.remove(e)
        logger.debug(f'chain_ends: {chain_ends}')
        flags=[]
        for p in combinations(chain_ends,2):
            n,m=p
            logger.debug(f'query chain end pair {n} {m}')
            if len(list(nx.shortest_path(G,n,m)))<7:
                logger.debug(f'removing {m} from chain_ends')
                flags.append(m)
        logger.debug(f'flags {flags}')
        for m in flags:
            chain_ends.remove(m)
        logger.debug(f'chain_ends: {chain_ends}')
        pathlength={}
        # if 0<len(chain_ends)<3:
        for e in ends:
            pathlength[e]={}
            for c in chain_ends:
                pathlength[e][c]=len(list(nx.shortest_path(G,e,c)))
        logger.debug(f'pathlength {pathlength}')
        avl=[]
        en=[]
        for k,v in pathlength.items():
            av=0.0
            if len(v)>0:
                for kk,vv in v.items():
                    av+=vv
                if len(v)>0:
                    avl.append(av/len(v))
                en.append(k)
        if avl:
            avl=np.array(avl)
            heads=[en[np.argmax(avl)]]
        
        tails=chain_ends
        shortest_paths=pathlength
        logger.debug(f'heads {heads} tails {tails} shortest_paths {shortest_paths}')
    return heads,tails,shortest_paths

def do_psfgen(resid,DB,lenfac=1.2,minimize_steps=500,sample_steps=5000,nsamples=10,sample_temperature=300):
    charmm_topfile=DB.get_charmm_topfile(resid)
    if charmm_topfile is None:
        return -2
    topo=DB.get_topo(resid)
    my_logger(f'{os.path.basename(topo.charmmtopfile)}',logger.info,just='^',frame='*',fill='*')
    heads,tails,shortest_paths=head_tail_atoms(topo)
    logger.debug(f'resi {resid} heads {heads} tails {tails} shortest_paths {shortest_paths}')
    tasklist=[
            {'restart':{'psf': f'{resid}-init.psf','pdb': f'{resid}-init.pdb'}},
            {'md': {'ensemble':'minimize','nsteps':0,'minimize':minimize_steps,'dcdfreq':0,'xstfreq':0,'temperature':100}},
            ]
    if shortest_paths and any([len(v)>0 for v in shortest_paths.values()]) and heads and tails and len(tails)<3:
        base_md={'ensemble':'NVT','nsteps':10000,'dcdfreq':100,'xstfreq':100,'temperature':100}
        groups={}
        if len(tails)<4:
            groups={'repeller':{'atomnames': [heads[0]]}}
        distances={}
        harmonics={}
        for i in range(len(tails)):
            groups[f'tail{i+1}']={'atomnames': [tails[i]]}
            if len(tails)<4:
                distances[f'repeller_tail{i+1}']={'groups': ['repeller',f'tail{i+1}']}
        dists={}
        for i in range(len(tails)):
            dists[i]=shortest_paths[heads[0]][tails[i]]*lenfac
        for i in range(len(tails)):
            for j in range(i+1,len(tails)):
                distances[f'tail{i+1}_tail{j+1}']={'groups': [f'tail{i+1}',f'tail{j+1}']}
                harmonics[f'tail{i+1}_tail{j+1}_attract']={
                        'colvars': [f'tail{i+1}_tail{j+1}'],
                        'forceConstant': 10.0,
                        'distance':[3.0]}
        if len(tails)<4:
            name='repeller_tail'+''.join(f'{i+1}' for i in range(len(tails)))
            colvars=[]               
            distance=[]
            for i in range(len(tails)):
                colvars.append(f'repeller_tail{i+1}')
                distance.append(dists[i])
            harmonics[name]={'colvars':colvars,'distance':distance,'forceConstant':10.0}
        colvar_specs={'groups':groups,'distances':distances,'harmonics':harmonics}
        base_md['colvar_specs']=colvar_specs
        assert 'minimize' not in base_md
        logger.debug(f'base_md {base_md}')
        tasklist.append({'md':base_md})
    else:
        logger.debug(f'**** Could not analyze structure of resi {resid}')
    if len(heads)>0:
        tasklist.append({'manipulate':{'mods':{'orient':[f'z,{heads[0]}']}}})
    tasklist.append({'md':{'ensemble':'NVT','nsteps':sample_steps,'dcdfreq':sample_steps//nsamples,'xstfreq':100,'temperature':sample_temperature,'index':99}})

    C=Controller(userspecs={'tasks':tasklist},quiet=True)
    # First we de-novo generate a pdb/psf file using seeded internal coordinates
    W=Psfgen(C.config)
    if not charmm_topfile in W.charmmff_config['standard']['topologies']:
        W.charmmff_config['standard']['topologies'].append(charmm_topfile)
    if charmm_topfile.endswith('detergent.str'):
        needed=[x for x in DB.all_charmm_topology_files if x.endswith('sphingo.str')][0]
        W.charmmff_config['standard']['topologies'].append(needed)
    if charmm_topfile.endswith('initosol.str'):
        needed=os.path.join(DB.toppardir,'stream/carb/toppar_all36_carb_glycolipid.str')
        W.charmmff_config['standard']['topologies'].append(needed)

    W.newscript('init')
    if topo.to_psfgen(W)==-1:
        logger.info(f'skipping {resid} -- no ICs')
        with open(f'{resid}-topo.rtf','w') as f:
            topo.to_file(f)
        return -1
    W.writescript(f'{resid}-init')
    W.runscript()
    with open('init.log','r') as f:
        loglines=f.read().split('\n')
    for l in loglines:
        if 'Warning: failed to guess coordinates' in l:
            logger.warning(f'Could not psfgen-build {resid}')
            return -2

    # now run the tasks to minimize, stretch, orient along z, equilibrate, and sample
    tasks=C.tasks
    for task in tasks:
        if task.taskname=='md':
            if charmm_topfile.endswith('str') and not charmm_topfile in task.writers['namd2'].charmmff_config['standard']['parameters']:
                task.writers['namd2'].charmmff_config['standard']['parameters'].append(charmm_topfile)
            if charmm_topfile.endswith('sphingo.str'):
                needed=[x for x in DB.all_charmm_topology_files if x.endswith('lps.str')][0]
                task.writers['namd2'].charmmff_config['standard']['parameters'].append(needed)
            if charmm_topfile.endswith('detergent.str'):
                needed=[x for x in DB.all_charmm_topology_files if x.endswith('sphingo.str')][0]
                task.writers['namd2'].charmmff_config['standard']['parameters'].append(needed)
                needed=[x for x in DB.all_charmm_topology_files if x.endswith('cholesterol.str')][0]
                task.writers['namd2'].charmmff_config['standard']['parameters'].append(needed)
            if charmm_topfile.endswith('inositol.str'):
                needed=os.path.join(DB.toppardir,'stream/carb/toppar_all36_carb_glycolipid.str')
                task.writers['namd2'].charmmff_config['standard']['parameters'].append(needed)

    result=C.do_tasks()
    for k,v in result.items():
        logger.debug(f'{k}: {v}')
        if v['result']!=0:
            return -1

    if os.path.exists('99-00-md-NVT.namd'):
        with open('99-00-md-NVT.namd','r') as f:
            lines=f.read().split('\n')
        par=[]
        for l in lines:
            if l.startswith('parameters'):
                if 'toppar/' in l:
                    par.append(l.split()[1].split('toppar/')[1])
        with open('parameters.txt','w') as f:
            for p in par:
                f.write(f'{p}\n')
        
    # now sample
    W=VMD(C.config)
    W.newscript('sample')
    W.addline(f'mol new {resid}-init.psf')
    W.addline(f'mol addfile 99-00-md-NVT.dcd waitfor all')
    W.addline(f'set a [atomselect top all]')
    W.addline(f'set ref [atomselect top all]')
    W.addline(f'$ref frame 0')
    W.addline(f'$ref move [vecscale -1 [measure center $ref]]')
    W.addline(r'for { set f 0 } { $f < [molinfo top get numframes] } { incr f } {')
    W.addline( '    $a frame $f')
    W.addline( '    $a move [measure fit $a $ref]')
    W.addline(f'    $a writepdb {resid}-[format %02d $f].pdb')
    W.addline(r'}')
    W.writescript()
    W.runscript()

    pdbs=[os.path.basename(x) for x in glob.glob(f'{resid}-??.pdb')]
    smin=[]
    smax=[]
    zmin=[]
    zmax=[]
    for pdb in pdbs:
        p=PDBParser(PDBcode=os.path.splitext(pdb)[0]).parse()
        z=np.array([x.z for x in p.parsed['ATOM']])
        s=np.array([x.serial for x in p.parsed['ATOM']])
        zmin_idx=np.argmin(z)
        zmax_idx=np.argmax(z)
        smin.append(s[zmin_idx])
        zmin.append(z[zmin_idx])
        smax.append(s[zmax_idx])
        zmax.append(z[zmax_idx])
    measures=pd.DataFrame({'pdb':pdbs,'bottom_serial':smin,'bottom_z':zmin,'top_serial':smax,'top_z':zmax})
    measures.to_csv('orientations.dat',header=True,index=False)
    return 0

def do_cleanup(resname,dirname):
    cwd=os.getcwd()
    os.chdir(dirname)
    files=glob.glob('*')
    files.remove('init.tcl')
    files.remove('orientations.dat')
    files.remove('parameters.txt')
    files.remove(f'{resname}-init.psf')
    for f in glob.glob(f'{resname}-*.pdb'):
        files.remove(f)
    for f in files:
        os.remove(f)
    os.chdir(cwd)

def do_resi(resi,DB,outdir='data',faildir='fails',force=False,lenfac=1.2,cleanup=True,minimize_steps=500,sample_steps=5000,nsamples=10,sample_temperature=300):
    cwd=os.getcwd()
    if (not os.path.exists(os.path.join(outdir,resi))) or force:
        os.mkdir('tmp')
        os.chdir('tmp')
        result=do_psfgen(resi,DB,lenfac=lenfac,minimize_steps=minimize_steps,sample_steps=sample_steps,nsamples=nsamples,sample_temperature=sample_temperature)
        os.chdir(cwd)
        if result==0:
            if cleanup: do_cleanup(resi,'tmp')
            shutil.move('tmp',os.path.join(outdir,resi))
        elif result==-2:
            my_logger(f'RESI {resi} is not found',logger.warning,just='^',frame='*',fill='*') 
        else:
            if os.path.exists(os.path.join(faildir,resi)):
                shutil.rmtree(os.path.join(faildir,resi))
            shutil.move('tmp',os.path.join(faildir,resi))
    else:
        logger.info(f'RESI {resi} in {outdir}; use \'--force\' to recalculate')

def make_RESI_database(args):
    streams=args.streams
    loglevel_numeric=getattr(logging,args.diagnostic_log_level.upper())
    if args.diagnostic_log_file:
        if os.path.exists(args.diagnostic_log_file):
            shutil.copyfile(args.diagnostic_log_file,args.diagnostic_log_file+'.bak')
        logging.basicConfig(filename=args.diagnostic_log_file,filemode='w',format='%(asctime)s %(name)s %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    DB=CharmmResiDatabase(streams=streams)
    outdir=args.output_dir
    faildir=args.fail_dir
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not os.path.exists(faildir):
        os.mkdir(faildir)
    if os.path.exists('tmp'):
        shutil.rmtree('tmp')
    resi=args.resi
    if resi is not None:
        my_logger(f'RESI {resi}',logger.info,just='^',frame='*',fill='*')
        do_resi(resi,DB,outdir=outdir,faildir=faildir,force=args.force,cleanup=args.cleanup,lenfac=args.lenfac,minimize_steps=args.minimize_steps,sample_steps=args.sample_steps,nsamples=args.nsamples,sample_temperature=args.sample_temperature)
    else:
        nresi=len(DB.charmm_resnames)
        for i,resi in enumerate(DB.charmm_resnames):
            my_logger(f'RESI {resi} ({i+1}/{nresi})',logger.info,just='^',frame='*',fill='*')
            do_resi(resi,DB,outdir=outdir,faildir=faildir,force=args.force,cleanup=args.cleanup,lenfac=args.lenfac,minimize_steps=args.minimize_steps,sample_steps=args.sample_steps,nsamples=args.nsamples,sample_temperature=args.sample_temperature)

