# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import networkx as nx
from collections import UserDict, UserList
import numpy as np
from .config import ResourceManager
from .controller import Controller
from .stringthings import linesplit
from .scriptwriters import Psfgen,VMD
import glob 
import os
import logging
import shutil

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
                return a.element
        # else:
        #     logger.debug(f'Could not find atom with name {a}')
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
        data=np.array(list(map(float,toks[6:])))
        if np.any(data==0.0):
            logger.debug(f'{toks[1:5]}: no ic data')
            logger.debug(f'{data}')
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
        datacards=lines[didx:]
        # logger.debug(f'{len(datacards)} datacards in resi {self.resname}')
        self.atomcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('ATOM') or datacards[i].startswith('Atom') ]
        self.groupcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('GROUP') or datacards[i].startswith('Group') ]
        self.bondcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('BOND') or datacards[i].startswith('Bond') ]
        self.doublecard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('DOUBLE') or datacards[i].startswith('Double') ]
        self.triplecard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('TRIPLE') or datacards[i].startswith('Triple') ]
        self.ICcard_idx=[i for i in range(len(datacards)) if datacards[i].startswith('IC') or datacards[i].startswith('Ic') ]
        self.atoms=CharmmTopAtomList([])
        self.atoms_in_group={i:CharmmTopAtomList([]) for i in range(len(self.groupcard_idx))}
        self.comment_strings=[]
        if self.groupcard_idx:
            for gn,gi in enumerate(self.groupcard_idx):
                dummy,comment=linesplit(datacards[gi])
                if len(comment.strip())>0:
                    self.comment_strings.append(comment)
                ai=1
                while datacards[gi+ai].startswith('ATOM') or datacards[gi+ai].startswith('Atom') :
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
            assert element1!='?',f'no element for atom {node1} in {self.resname} {self.charmmtopfile}'
            assert element2!='?',f'no element for atom {node2} in {self.resname} {self.charmmtopfile}'
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
    
    def to_psfgen(self,W):
        for ic in self.IC:
            if not ic.empty:
                if ic.atoms[1]==self.atoms[0].name and ic.dihedral_type != 'improper':
                    break
        else:
            logger.debug(f'atom {self.atoms[0].name} is not represented in the ic list')
            return -1
        W.addline(r'segment A {')
        W.addline(r'    first none')
        W.addline(r'    last none')
        W.addline(f'    resid 1 {self.resname}')
        W.addline(r'}')
        W.addline(f'psfset coord A 1 {self.atoms[0].name} '+r'{0 0 0}')
        b1a1,b1a2,b1l=ic.bonds[0].name1,ic.bonds[0].name2,ic.bonds[0].length
        partner=b1a1
        assert b1a1!=self.atoms[0].name
        W.addline(f'psfset coord A 1 {partner} '+r'{'+f'{b1l:.5f}'+r' 0 0}')
        a1a1,a1a2,a1a3,a1d=ic.angles[0].name1,ic.angles[0].name2,ic.angles[0].name3,ic.angles[0].degrees
        assert a1a1==partner
        assert a1a2==self.atoms[0].name,f'{self.atoms[0].name} {ic.atoms}: {ic.angles[0].name1} {ic.angles[0].name2} {ic.angles[0].name3}'
        assert a1a3!=partner
        W.addline(f'set a [expr {a1d}*acos(-1.0)/180.0]')
        W.addline(r'set x [expr 1.5283*cos($a)]')
        W.addline(r'set y [expr 1.5283*sin($a)]')
        W.addline(f'psfset coord A 1 {a1a3} [list $x $y 0.0]')
        return 0

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
        lines=f.read().upper().split('\n')
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
        toppardir=self.resources['charmmff']['toppar']
        rtftops=glob.glob(os.path.join(toppardir,'top_*.rtf'))
        M=[]
        for f in rtftops:
            M.extend(getMasses(f))
        strtops=[]
        tops=[]
        for stream in streams:
            tops+=glob.glob(os.path.join(toppardir,f'*{stream}*.rtf'))
            streamdir=os.path.join(toppardir,f'stream/{stream}')
            strfiles=glob.glob(os.path.join(streamdir,'*.str'))
            idx=['list' in x for x in strfiles].index(True)
            strfiles.remove(strfiles[idx])
            strtops+=glob.glob(os.path.join(streamdir,'*.str'))
        for f in strtops:
            M.extend(getMasses(f))
        logger.debug(f'{len(M)} mass records')
        self.M=CharmmMasses(M)

        data={}
        for t in tops+strtops:
            sublist=getResis(t,self.M)
            for resi in sublist:
                resi.charmmtopfile=t
                data[resi.resname]=resi
        super().__init__(data)
        logger.debug(f'{len(self)} CHARMM topology residues parsed')
            # logger.debug(f'{os.path.basename(t)} {len(sublist)} RESIs')
        self.charmm_resnames=list(sorted(list(self.keys())))
    
    def get_charmm_topfile(self,charmm_resid):
        elem=self.get(charmm_resid,None)
        if elem is not None:
            return elem.charmmtopfile
        return ''
    
    def get_topo(self,charmm_resid):
        elem=self.get(charmm_resid,None)
        if elem is not None:
            return elem
        return None

def do_psfgen(resid,DB,cleanup=True):
    charmm_topfile=DB.get_charmm_topfile(resid)
    topo=DB.get_topo(resid)
    i=0
    toprefatom=topo.atoms[i]
    while toprefatom.element=='H':
        i+=1
        toprefatom=topo.atoms[i]
    toprefatom=toprefatom.name
    logger.debug(f'{resid} topatom {toprefatom}')
    g=topo.to_graph(includeH=False)
    ends=[]
    mid=[]
    for n in g.__iter__():
        nn=len(list(g.neighbors(n)))
        if nn==1:
            ends.append(n)
        elif nn==3:
            mid.append(n)
    # logger.debug(f'{resid} ends {ends} mid {mid}')
    eb=[]
    for e in ends:
        L=nx.shortest_path(g,e,toprefatom)
        eb.append(len(L))
    eb=np.array(eb)
    idx=eb.argsort()
    eb_sorted=eb[idx]
    idx2=eb_sorted>20
    ends_sorted=np.array(ends)[idx]
    top_ends=ends_sorted[idx2]
    logger.debug(f'top ends {top_ends}')

    tasklist=[{'restart':{
                        'psf': f'{resid}-init.psf',
                        'pdb': f'{resid}-init.pdb'
                        },
                    },    
                    {'md': {'ensemble':'minimize','nsteps':0,'minimize':1000,'dcdfreq':0,'xstfreq':0,'temperature':100}}]
    if len(top_ends)==2:
        top_dist=eb_sorted[idx2].mean()*30.0/25.0
        tasklist.append({'md':{'ensemble':'NVT','nsteps':10000,'dcdfreq':100,'xstfreq':100,'temperature':100,
                        'colvar_specs': {
                            'groups': {
                                'repeller': {'atomnames': [toprefatom]},
                                'tail1': {'atomnames': [top_ends[0]]},
                                'tail2': {'atomnames': [top_ends[1]]}
                                },
                            'distances': {
                                'tail1_tail2': {'groups': ['tail1','tail2']},
                                'repeller_tail1': {'groups': ['repeller','tail1']},
                                'repeller_tail2':{'groups': ['repeller','tail2']}
                                },
                            'harmonics': {
                                'tail1_tail2_attract': {
                                    'colvars': ['tail1_tail2'],
                                    'forceConstant': 10.0,
                                    'distance':[3.0]},
                                'repeller_tail12_repel': {
                                    'colvars': ['repeller_tail1','repeller_tail2'],
                                    'forceConstant': 10.0,
                                    'distance':[top_dist,top_dist]}
                            }
                        }}})
    tasklist.append({'manipulate':{'mods':{'orient':[f'z,{toprefatom}']}}})
    tasklist.append({'md':{'ensemble':'NVT','nsteps':5000,'dcdfreq':500,'xstfreq':100,'temperature':300,'index':99}})

    C=Controller(userspecs={'tasks':tasklist})
    # First we de-novo generate a pdb/psf file using seeded internal coordinates
    W=Psfgen(C.config)
    if not charmm_topfile in W.charmmff_config['standard']['topologies']:
        W.charmmff_config['standard']['topologies'].append(charmm_topfile)
    W.newscript('init')
    if topo.to_psfgen(W)==-1:
        logger.info(f'skipping {resid} -- no ICs')
        return -1
    W.writescript(f'{resid}-init')
    W.runscript()

    # now run the tasks to minimize, stretch, orient along z, equilibrate, and sample
    tasks=C.tasks
    for task in tasks:
        if task.taskname=='md':
            if charmm_topfile.endswith('str') and not charmm_topfile in task.writers['namd2'].charmmff_config['standard']['parameters']:
                task.writers['namd2'].charmmff_config['standard']['parameters'].append(charmm_topfile)
    C.do_tasks()

    # now sample and clean up
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
    W.addline(f'    $a writepdb {resid}-$f.pdb')
    W.addline(r'}')
    W.writescript()
    W.runscript()

    if cleanup:
        files=glob.glob('*')
        files.remove('init.tcl')
        files.remove(f'{resid}-init.psf')
        for f in glob.glob(f'{resid}-*.pdb'):
            files.remove(f)
        for f in glob.glob('*manipulate*pdb'):
            files.remove(f)
        for f in files:
            os.remove(f)

    return 0

def do_resi(resi,DB):
    cwd=os.getcwd()
    if not os.path.exists(os.path.join('data',resi)):
        os.mkdir('tmp')
        os.chdir('tmp')
        result=do_psfgen(resi,DB)
        os.chdir(cwd)
        if result==0:
            shutil.move('tmp',os.path.join('data',resi))
        else:
            shutil.rmtree('tmp')
    else:
        logger.info(f'...done')

def make_pdb_database(args):
    streams=args.streams
    loglevel_numeric=getattr(logging,args.loglevel.upper())
    logging.basicConfig(format='%(levelname)s> %(message)s',level=loglevel_numeric)
    DB=CharmmResiDatabase(streams=streams)
    if os.path.exists('tmp'):
        shutil.rmtree('tmp')
    resi=args.resi
    if resi is not None:
        do_resi(resi,DB)
    else:
        for resi in DB.charmm_resnames:
            do_resi(resi,DB)

