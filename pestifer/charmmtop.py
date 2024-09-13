# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Facilitates reading and parsing of CHARMM RESI and MASS records from rtp and str files
#
import glob 
import logging
import os
import networkx as nx
import numpy as np
from collections import UserDict, UserList
from itertools import combinations, compress, batched
from .config import ResourceManager
from .stringthings import linesplit

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

def charmmBonds(card:str,degree=1):
    result=CharmmBondList([])
    toks=[x.strip() for x in card.split()][1:]
    for p in batched(toks,2):
        result.append(CharmmBond(p,degree=degree))
    return result

class CharmmBond:
    def __init__(self,names,degree=1):
        self.name1,self.name2=names
        self.degree=degree

    def __eq__(self,other):
        c1=self.name1==other.name1 and self.name2==other.name2
        c2=self.name1==other.name2 and self.name2==other.name1
        return c1 or c2
    
class CharmmBondList(UserList):
    def __init__(self,data):
        self.data=data

def charmmAngles(card:str):
    result=CharmmAngleList([])
    toks=[x.strip() for x in card.split()][1:]
    for p in batched(toks,3):
        result.append(CharmmAngle(p))
    return result

class CharmmAngle:
    def __init__(self,names):
        self.name1,self.name2,self.name3=names

    def __eq__(self,other):
        r1=self.name2==other.name2
        c1=self.name1==other.name1 and self.name3==other.name3
        c2=self.name1==other.name3 and self.name3==other.name1
        return r1 and (c1 or c2)

class CharmmAngleList(UserList):
    def __init__(self,data):
        self.data=data

def charmmDihedrals(card,dihedral_type='proper'):
    result=CharmmDihedralList([])
    toks=[x.strip() for x in card.split()][1:]
    for p in batched(toks,4):
        result.append(CharmmDihedral(p,dihedral_type=dihedral_type))
    return result

class CharmmDihedral:
    def __init__(self,names,dihedral_type='proper'):
        self.name1,self.name2,self.name3,self.name4=names
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
        toks=[x.strip() for x in ICstring.split()]
        data=np.array(list(map(float,toks[5:])))
        if data[0]==0.0 or data[1]==0.0 or data[3]==0.0 or data[4]==0.0:
            # logger.debug(f'{toks[1:5]}: missing ic data: {data}')
            self.empty=True
            # return
        self.atoms=toks[1:5]
        if self.atoms[2].startswith('*'):
            self.atoms[2]=self.atoms[2][1:]
            self.dihedral_type='improper'
        else:
            self.dihedral_type='proper'
        if self.dihedral_type=='proper':
            b0=CharmmBond([self.atoms[0],self.atoms[1]])
            b0.length=float(toks[5])
            b1=CharmmBond([self.atoms[2],self.atoms[3]])
            b1.length=float(toks[9])
            self.bonds=CharmmBondList([b0,b1])
            a1=CharmmAngle([self.atoms[0],self.atoms[1],self.atoms[2]])
            a1.degrees=float(toks[6])
            a2=CharmmAngle([self.atoms[1],self.atoms[2],self.atoms[3]])
            a2.degrees=float(toks[8])
            self.angles=CharmmAngleList([a1,a2])
            self.dihedral=CharmmDihedral([self.atoms[0],self.atoms[1],self.atoms[2],self.atoms[3]])
            self.dihedral.degrees=float(toks[7])
        else:
            b0=CharmmBond([self.atoms[0],self.atoms[2]])
            b0.length=float(toks[5])
            b1=CharmmBond([self.atoms[2],self.atoms[3]])
            b1.length=float(toks[9])
            self.bonds=CharmmBondList([b0,b1])
            a1=CharmmAngle([self.atoms[0],self.atoms[2],self.atoms[3]])
            a1.degrees=float(toks[6])
            a2=CharmmAngle([self.atoms[1],self.atoms[2],self.atoms[3]])
            a2.degrees=float(toks[8])
            self.angles=CharmmAngleList([a1,a2])
            self.dihedral=CharmmDihedral([self.atoms[0],self.atoms[1],self.atoms[2],self.atoms[3]],dihedral_type='improper')
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
        self.error_code=0
        lines=[x.strip() for x in blockstring.split('\n')]
        titlecard=lines[0]
        assert titlecard.startswith('RESI'),f'bad title card in RESI: [{titlecard}]'
        titledata,titlecomment=linesplit(titlecard)
        tctokens=titledata.split()
        self.resname=tctokens[1]
        self.charge=float(tctokens[2])
        self.synonym=titlecomment.strip()

        didx=1
        while lines[didx].startswith('!'): didx+=1
        cards=[x.upper() for x in lines[didx:]]
        # logger.debug(f'first post-title card {cards[0]}')
        comments=[]
        datacards=[]
        for c in cards:
            rawcarddata,cardcomment=linesplit(c)
            carddata=rawcarddata.lstrip().upper()  # no indent, no case-sensitivity
            if len(carddata)>0:    datacards.append(carddata)
            if len(cardcomment)>0: comments.append(cardcomment)
        # logger.debug(f'{self.resname}: {len(datacards)} datacards and {len(comments)} comments')
        isatomgroup=[d.startswith('ATOM') or d.startswith('GROU') for d in datacards]
        atomgroupcards=compress(datacards,isatomgroup)
        self.atoms=CharmmTopAtomList([])
        self.atoms_in_group={}
        g=0
        for card in atomgroupcards:
            if card.startswith('ATOM'):
                a=CharmmTopAtom(card,masses=masses)
                self.atoms.append(a)
                if not g in self.atoms_in_group:
                    self.atoms_in_group[g]=CharmmTopAtomList([])
                self.atoms_in_group[g].append(a)
            elif card.startswith('GROU'):
                g+=1
        # logger.debug(f'{len(self.atoms)} atoms processed in {len(self.atoms_in_group)} groups')
        self.atomdict={a.name:a for a in self.atoms}
        isbond=[d.startswith('BOND') for d in datacards]
        bondcards=compress(datacards,isbond)
        self.bonds=CharmmBondList([])
        for card in bondcards:
            self.bonds.extend(charmmBonds(card))
        bondcards=None
        isdouble=[d.startswith('DOUB') for d in datacards]
        doublecards=compress(datacards,isdouble)
        for card in doublecards:
            self.bonds.extend(charmmBonds(card,degree=2))
        istriple=[d.startswith('TRIP') for d in datacards]
        triplecards=compress(datacards,istriple)
        for card in triplecards:
            self.bonds.extend(charmmBonds(card,degree=3))
        isIC=[d.startswith('IC') for d in datacards]
        ICcards=compress(datacards,isIC)
        self.IC=[]
        ic_atom_names=[]
        for card in ICcards:
            IC=CharmmTopIC(card)
            if any([x not in self.atomdict.keys() for x in IC.atoms]):
                self.error_code=-5
            else:
                ic_atom_names.extend(IC.atoms)
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
            # logger.debug(f'bond {node1} {node2}')
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
    
    def head_tail_atoms(self):
        """ Identify head and tail atoms """
        G=self.to_graph(includeH=False)
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
            # logger.debug(f'edge_partner_index {edge_partner_idx}')
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
                # head atom is the hydroxyl O on ringA
                for n in nx.neighbors(G,a):
                    if G.nodes[n]['element']=='O':
                        heads=[n]
            paths=[]
            for a in G.__iter__():
                paths.append(len(nx.shortest_path(G,a,heads[0])))
            # tail is atom furthest away from head
            paths=np.array(paths)
            l_idx=np.argsort(paths)[-1]
            tails=[list(G.__iter__())[l_idx]]
        else:
            # all atoms with only a single neighbor are "ends"' remember there are no H's here
            ends=[]
            for n in G.__iter__():
                nn=len(list(G.neighbors(n)))
                if nn==1:
                    ends.append(n)
            # for each end-neighbor, determine if any of them are bound to the same neighbor
            ecount={}
            for e in ends:
                nn=list(G.neighbors(e))[0]
                if not nn in ecount:
                    ecount[nn]=[]
                ecount[nn].append(e)
            logger.debug(f'ecount {ecount}')
            # for any atom that is a common neighbor to two or more "ends", it becomes a new
            # end and the old ends are removed from the "ends"
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
    
    def to_file(self,f):
        f.write(self.blockstring)

    def show_error(self,buf):
        if self.error_code==-5:
            buf(f'{self.resname} references atoms in IC\'s that are not declared at ATOM\'s')

    def to_psfgen(self,W,**kwargs):
        refic_idx=kwargs.get('refic-idx',0)
        # find the first heavy atom that is atom-i in an IC w/o H that is not an improper
        refatom=None
        refic=None
        logger.debug(f'trying to find a reference atom from {len(self.atoms)} atoms in this resi')
        i=0
        try:
            is_proper=[((not x.empty) and (x.dihedral_type=='proper') and (not any([self.atomdict[a].element=='H' for a in x.atoms]))) for x in self.IC]
        except:
            logger.warning(f'Could not parse IC\'s')
            return -1
        workingIC=list(compress(self.IC,is_proper))
        for ic in workingIC:
            logger.debug(f'{str(ic)}')
        if not workingIC:
            logger.warning(f'No valid IC for {self.resname}')
            return -1
        nWorkingIC=len(workingIC)
        logger.debug(f'{nWorkingIC}/{len(self.IC)} IC\'s with proper dihedrals and no hydrogens')
        refic=workingIC[refic_idx]
        refatom=self.atomdict[refic.atoms[1]]
        if refatom is None:
            logger.debug(f'Unable to identify a reference atom')
            return -1
        logger.debug(f'refIC is {str(refic)}')
        logger.debug(f"{[self.atomdict[a].element=='H' for a in refic.atoms]}")
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
        # assert b1a1!=refatom.name
        # assert b1a2==refatom.name
        # partner is the first atom in the IC, put along x-axis
        W.addline(f'psfset coord A 1 {partner} '+r'{'+f'{b1l:.5f}'+r' 0 0}')
        a1a1,a1a2,a1a3,a1d=refic.angles[0].name1,refic.angles[0].name2,refic.angles[0].name3,refic.angles[0].degrees
        # assert a1a1==partner
        # assert a1a2==refatom.name,f'{refatom.name} {refic.atoms}: {refic.angles[0].name1} {refic.angles[0].name2} {refic.angles[0].name3}'
        # assert a1a3!=partner
        W.addline(f'set a [expr {a1d}*acos(-1.0)/180.0]')
        W.addline(f'set x [expr {b1l}*cos($a)]')
        W.addline(f'set y [expr {b1l}*sin($a)]')
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
        # logger.debug(f'resi {resi.resname} parsed')
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
    def __init__(self):  # initialize from standard topology files in the topmost directory
        self.resources=ResourceManager()
        self.toppardir=self.resources['charmmff']['toppar']
        self.customdir=self.resources['charmmff']['custom']
        self.toppath=[self.toppardir,self.customdir]
        # get all atom types from standard charmmff topology files
        stdtops=glob.glob(os.path.join(self.toppardir,'top_*.rtf'))
        stdtops=[x for x in stdtops if 'ljpme' not in x]
        stdtops.extend(glob.glob(os.path.join(self.toppardir,'toppar*.str')))
        self.all_charmm_topology_files=stdtops
        M=[]
        for f in stdtops:
            M.extend(getMasses(f))
        
        self.M=CharmmMasses(M)
        self.charmm_resnames=[]
        data={}
        for f in stdtops:
            sublist=getResis(f,self.M)
            for resi in sublist:
                resi.charmmtopfile=f
                if f.endswith('.rtf'): # "top_allXX_<streamname>.rtf"
                    stream=os.path.splitext(os.path.basename(f))[0].split('_')[2]
                elif f.endswith('.str'): # a stream file in the base directory, most likely toppar_water_ions.str
                    stream='_'.join(os.path.splitext(os.path.basename(f))[0].split('_')[1:])
                if not stream in data:
                    data[stream]={}
                data[stream][resi.resname]=resi
                self.charmm_resnames.append(resi.resname)

        super().__init__(data)
        logger.debug(f'{len(stdtops)} CHARMM topology files scanned')
        logger.debug(f'Streams initiated: {" ".join(list(data.keys()))}')
        self.charmm_resnames.sort()
        logger.debug(f'{len(self.charmm_resnames)} RESI\'s parsed')
        logger.debug(f'{len(self.M)} atomic mass records parsed')

    def add_topology(self,topfile,streamnameoverride=''):
        if not topfile.startswith('/'):
            for p in self.toppath:
                top_abs=os.path.join(p,topfile)
                if os.path.exists(top_abs):
                    break
            else:
                logger.debug(f'{topfile} not found in {self.toppath}')
                return -1
        else:
            top_abs=topfile
        if top_abs in self.all_charmm_topology_files:
            logger.debug(f'{topfile} has already been incorporated')
            return self

        self.all_charmm_topology_files.append(topfile)
        newM=CharmmMasses(getMasses(topfile))
        self.M.update(newM)
        sublist=getResis(topfile,self.M)
        for resi in sublist:
            resi.charmmtopfile=topfile
            if topfile.endswith('.rtf') or (topfile.startswith('stream/') and topfile.endswith('.str')):
                stream=os.path.splitext(os.path.basename(topfile))[0].split('_')[2] # "toppar_allXX_<streamname>[_YYY].[str,rtf]"
            else:
                stream='_'.join(os.path.splitext(os.path.basename(topfile))[0].split('_')[1:])
            if streamnameoverride != '':
                stream=streamnameoverride
            if not stream in self:
                self[stream]={}
            self[stream][resi.resname]=resi
            self.charmm_resnames.append(resi.resname)

        self.charmm_resnames.sort()
        logger.debug(f'Appended {len(sublist)} RESI\'s and {len(newM)} mass records')
        logger.debug(f'Current streams: {" ".join(list(self.keys()))}')
        return self

    def add_stream(self,streamname):
        streamdir=os.path.join(self.toppardir,f'stream/{streamname}')
        if not os.path.exists(streamdir):
            logger.debug(f'{streamdir} not found')
            return self

        allstream_strs=glob.glob(os.path.join(streamdir,'*.str'))
        include_strs=[not ('cationpi_wyf' in x or 'list' in x or 'model' in x or 'lipid_prot' in x or 'ljpme' in x or 'llo' in x) for x in allstream_strs]
        stream_strs=list(compress(allstream_strs,include_strs))
        M=[]
        for f in stream_strs: # already got masses from all rtf files
            M.extend(getMasses(f))
        
        self.M.update(CharmmMasses(M))
        for t in stream_strs:
            sublist=getResis(t,self.M)
            for resi in sublist:
                if resi.resname in self.charmm_resnames:
                    logger.debug(f'RESI {resi.resname} is already in the database')
                else:
                    relt=t
                    if 'toppar/' in t:
                        relt=t.split('toppar/')[1]
                    resi.charmmtopfile=relt
                    if not streamname in self:
                        self[streamname]={}
                    self[streamname][resi.resname]=resi
                    self.charmm_resnames.append(resi.resname)
        return self

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
        logger.debug(f'resi {charmm_resid} is not found in any streams')
        return None


