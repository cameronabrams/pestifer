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
from itertools import compress, batched
from .resourcemanager import ResourceManager
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

    def get_mass(self,name):
        a=self.get_atom(name)
        if a is not None:
            if hasattr(a,'mass'):
                return a.mass
        return 0.0

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
        self.empty=False
        ICstring,dummy=linesplit(ICstring)
        toks=[x.strip() for x in ICstring.split()]
        data=np.array(list(map(float,toks[5:])))
        if data[0]==0.0 or data[1]==0.0 or data[3]==0.0 or data[4]==0.0:
            # logger.debug(f'{toks[1:5]}: missing ic data: {data}')
            self.empty=True
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

    def __str__(self):
        ans=','.join([f'{i+1}{self.atoms[i]}' for i in range(len(self.atoms))])
        ans+=':'
        ans+=','.join([f'{b.length:.4f}' for b in self.bonds])
        ans+=':'
        ans+=','.join([f'{a.degrees:.4f}' for a in self.angles])
        ans+=':'
        ans+=f'{self.dihedral.degrees}-{self.dihedral.dihedral_type}'
        if self.empty:
            ans+='-empty'
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
    
    def lipid_annotate(self):
        """ Identify head and tail atoms """
        m=self.metadata
        if m['stream']=='lipid' and m['substream']=='cholesterol':
            self.sterol_annotate()
        elif m['stream']=='lipid' and m['substream']=='detergent':
            self.detergent_annotate()
        else:
            self.generic_lipid_annotate()
    
    def sterol_annotate(self):
        self.annotation={}
        G=self.to_graph(includeH=False)
        heads=[]
        for a in G:
            # head atom is the carbon to which the hydroxyl O is bound
            for n in nx.neighbors(G,a):
                if G.nodes[n]['element']=='O':
                    heads=[a]
                    break
        paths=[]
        for a in G.__iter__():
            paths.append(len(nx.shortest_path(G,a,heads[0])))
        # tail is atom furthest away from head
        paths=np.array(paths)
        l_idx=np.argsort(paths)[-1]
        tails=[list(G.__iter__())[l_idx]]
        logger.debug(f'heads {heads} tails {tails}')    
        self.annotation['heads']=heads
        self.annotation['tails']=tails

    def detergent_annotate(self):
        self.annotation={}
        G=self.to_graph(includeH=False)
        # find the bond which, if deleted, produces two fragments such that
        # 1. one fragment contains only carbons (tail)
        # 2. the other fragment is as small as possible and contains at least one O (head)
        minheadsize=1000
        mindata={}
        for f in nx.bridges(G):
            logger.debug(f'f {f}')
            g=G.copy()
            g.remove_edge(*f)
            c=list(nx.connected_components(g))
            # if len(c)!=2: continue
            e=[[G.nodes[x]['element'] for x in y] for y in c]
            logger.debug(f'e {e}')
            allc=[all([x=='C' for x in y]) for y in e]
            haso=[any([x=='O' for x in y]) for y in e]
            logger.debug(f'allc {allc}')
            logger.debug(f'haso {haso}')
            if any(allc) and any(haso):
                tailidx=allc.index(True)
                headidx=1-tailidx
                logger.debug(f'tailidx {tailidx} headidx {headidx}')
                tailsize=len(e[tailidx])
                headsize=len(e[headidx])
                logger.debug(f'tailsize {tailsize} headsize {headsize}')
                if headsize<minheadsize:
                    minheadsize=headidx
                    mindata=dict(headg=list(c)[headidx],tailg=list(c)[tailidx])
        logger.debug(f'mindata {mindata}')
        headg=list(mindata['headg'])
        tailg=list(mindata['tailg'])
        # head reference is heaviest atom in head
        headmasses=np.array([self.atoms.get_mass(n) for n in headg])
        headatom=headg[np.argmax(headmasses)]
        logger.debug(f'headmasses {headmasses} headatom {headatom}')

        self.annotation['heads']=[headatom]

        maxpathlength=-1
        tailatom=None
        for ta in tailg:
            path=nx.shortest_path(G,source=headatom,target=ta)
            pathlength=len(path)
            if pathlength>maxpathlength:
                maxpathlength=pathlength
                tailatom=ta
        assert tailatom!=None
        self.annotation['tails']=[tailatom]

        logger.debug(f'annotation {self.annotation}')

    def generic_lipid_annotate(self):
        self.annotation={}
        G=self.to_graph(includeH=False)
        WG=G.copy()

        # find all carbon chains
        ctails=[]
        cbonds=[]
        hastails=True
        while hastails:
            hastails=False
            maxtailsize=-1
            result={}
            for f in nx.bridges(WG):
                g=WG.copy()
                g.remove_edge(*f)
                S=[g.subgraph(c).copy() for c in nx.connected_components(g)]
                c=list(nx.connected_components(g))
                if len(c)!=2: continue
                e=[[G.nodes[x]['element'] for x in y] for y in c]
                # logger.debug(f'e {e}')
                allc=[all([x=='C' for x in y]) for y in e]
                # haso=[any([x=='O' for x in y]) for y in e]
                # logger.debug(f'allc {allc}')
                # logger.debug(f'haso {haso}')
                if any(allc):
                    tail=c[allc.index(True)]
                    tailsize=len(tail)
                    # logger.debug(f'tail {tail} tailsize {tailsize}')
                    if tailsize>maxtailsize:
                        result=dict(tail=tail,nontailg=S[allc.index(False)],bond=f)
                        tailatoms=list(tail)
                        maxtailsize=tailsize
            if result:
                hastails=True
                # ignore methyls
                if len(tailatoms)>1:
                    ctails.append(tailatoms)
                    cbonds.append(result['bond'])
                    logger.debug(f'tailatoms {tailatoms}')
                WG=result['nontailg']

        logger.debug(f'ctails {ctails}')
        logger.debug(f'cbonds {cbonds}')
        WG=G.copy()

        # find the ends of all the carbon tails
        tails=[tail[[len(list(nx.neighbors(G,n))) for n in tail].index(1)] for tail in ctails]
        logger.debug(f'tailn {tails}')
        if len(tails)==2: # this is a "standard" lipid with two fatty acid tails
            queryheads=[k for k,n in WG.nodes.items() if n['element']!='C' and n['element']!='O']
            logger.debug(f'queryheads {queryheads}')
            if len(queryheads)==0: # only non-carbons in head are O's
                queryheads=[k for k,n in WG.nodes.items() if n['element']=='O']
            maxpathlen=[-1]*len(cbonds)
            claimedhead=['']*len(cbonds)
            for nc in queryheads:
                for i,b in enumerate(cbonds):
                    path=nx.shortest_path(WG,source=b[0],target=nc)
                    pathlen=len(path)
                    if pathlen>maxpathlen[i]:
                        maxpathlen[i]=pathlen
                        claimedhead[i]=nc
            logger.debug(f'claimedhead {claimedhead}')
            headcontenders=list(set(claimedhead))
            headvotes=np.array([claimedhead.count(x) for x in headcontenders])
            heads=[headcontenders[np.argmax(headvotes)]]
            logger.debug(f'heads {heads}')
        else:
            # this could be a wacky multiple-chain lipid with an extended head group
            # like cardiolipin or one of those bacterial glycolipids
            OG=WG.copy()
            for tail in ctails:
                for atom in tail:
                    OG.remove_node(atom)
            mins=len(OG)**2
            headatom=None
            for atom in OG:
                g=OG.copy()
                g.remove_node(atom)
                c=[len(x) for x in list(nx.connected_components(g))]
                while 1 in c:
                    c.remove(1)
                c=np.array(c)
                # logger.debug(f'testing atom {atom} -> {c}')
                if len(c)>1:
                    s=c.std()
                    if s<mins:
                        # logger.debug(f'c {c} {atom} {s}<{mins}')
                        headatom=atom
                        mins=s
            assert headatom!=None
            logger.debug(f'headatom {headatom}')
            heads=[headatom]

        self.annotation['heads']=heads
        self.annotation['tails']=tails
        self.annotation['shortest_paths']={head:{tail:len(nx.shortest_path(WG,source=head,target=tail)) for tail in tails} for head in heads}
    
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
        # partner is the first atom in the IC, put along x-axis
        W.addline(f'psfset coord A 1 {partner} '+r'{'+f'{b1l:.5f}'+r' 0 0}')
        a1a1,a1a2,a1a3,a1d=refic.angles[0].name1,refic.angles[0].name2,refic.angles[0].name3,refic.angles[0].degrees
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
    f,ext=os.path.splitext(os.path.basename(topfile))
    if not ext in ['.rtf','.str']:
        logger.debug(f'{topfile} extension \'{ext}\' not recognized')
        return None
    ftok=f.split('_')
    fcontent_type=ftok[0]
    if not fcontent_type in ['top','toppar']:
        logger.debug(f'{topfile} content type \'{fcontent_type}\' not recognized')
        return None
    stream,substream='',''
    idx=2 if 'all' in ftok[1] else 1
    if fcontent_type=='toppar':
        stream=ftok[idx]
        if len(ftok)==idx+1:
            substream=''
        else:
            substream=ftok[idx+1]
        if stream=='water' and substream=='ions':
            stream='water_ions'
            substream=''
    else:
        stream,substream=ftok[idx],''
    logger.debug(f'topfile {topfile} stream \'{stream}\' substream \'{substream}\'')
    R=[]
    metadat=dict(charmmtopfile=topfile,stream=stream,substream=substream)
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
            endidx=residx[-1]+i
            break
    bufs.append('\n'.join(lines[residx[-1]:endidx]))
    for block in bufs:
        resi=CharmmTopResi(block,masses=masses)
        resi.metadata=metadat
        R.append(resi)
    return R

class CharmmResiDatabase(UserDict):
    def __init__(self):  # initialize from standard topology files in the topmost directory
        self.resources=ResourceManager()
        self.toppardir=self.resources.get_charmmff_toppardir()
        self.customdir=self.resources.get_charmmff_customdir()
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
                stream,substream=resi.metadata['stream'],resi.metadata['substream']
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
            stream,substream=resi.metadata['stream'],resi.metadata['substream']
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
            if not t in self.all_charmm_topology_files:
                self.all_charmm_topology_files.append(t)
            sublist=getResis(t,self.M)
            for resi in sublist:
                if resi.resname in self.charmm_resnames:
                    logger.debug(f'RESI {resi.resname} is already in the database')
                else:
                    stream,substream=resi.metadata['stream'],resi.metadata['substream']
                    if not stream in self:
                        self[stream]={}
                    self[stream][resi.resname]=resi
                    self.charmm_resnames.append(resi.resname)
        return self

    def get_charmm_topfile(self,charmm_resid):
        if charmm_resid not in self.charmm_resnames:
            return None
        for stream,data in self.items():
            elem=data.get(charmm_resid,None)
            if elem is not None:
                return elem.metadata['charmmtopfile']
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


