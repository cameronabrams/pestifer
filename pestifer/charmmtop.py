# Author: Cameron F. Abrams, <cfa22@drexel.edu>
#
# Facilitates reading and parsing of CHARMM RESI and MASS records from rtp and str files
#
import logging
import networkx as nx
import numpy as np
from collections import UserDict, UserList
from itertools import compress, batched
from .stringthings import linesplit #, my_logger

logger=logging.getLogger(__name__)

class CharmmMassRecord:
    def __init__(self,record):
        tok=record
        tok,self.com=linesplit(record)
        tokens=[t.strip() for t in tok.split()]
        self.atom_type=tokens[2].upper()
        self.atom_mass=float(tokens[3])
        self.atom_element='?'
        try:
            self.atom_element=tokens[4]
        except:
            firstchar=self.atom_type[0]
            if not firstchar.isdigit():
                self.atom_element=self.atom_type[0].upper()
    
    def __str__(self):
        ans=f'{self.atom_type} {self.atom_mass:.4f}'
        if self.atom_element!='?':
            ans+=f' {self.atom_element}'
        if self.com:
            ans+=' '+self.com
        return ans

class CharmmMasses(UserDict):
    def __init__(self,mList):
        self.data={}
        for m in mList:
            self.data[m.atom_type]=m

class CharmmTopAtom:
    def __init__(self,atomstring):
        tokens=atomstring.split()
        self.name=tokens[1].upper()
        self.type=tokens[2].upper()
        self.charge=float(tokens[3])
        self.mass=0.0
        self.element='?'
        if self.name[0].isdigit():
            # logger.debug(f'Atom name {self.name} starts with a digit; setting inpatch')
            self.inpatch=True
    
    def __str__(self):
        ans=f'ATOM {self.name} {self.type} {self.charge:.4f}'
        if hasattr(self,'mass') and self.mass>0.0:
            ans+=f' {self.mass:.4f}'
        if hasattr(self,'element') and self.element!='?':
            ans+=f' {self.element}'
        return ans

    def __eq__(self,other):
        return self.name==other.name
    
    def set_mass(self,masses):
        """ Set the mass of this atom from a CharmmMasses object """
        # logger.debug(f'setting mass for {self.name} ({self.type}) using type {type(masses)}')
        if isinstance(masses,CharmmMasses):
            m=masses.get(self.type,None)
            if m is not None:
                self.mass=m.atom_mass
                self.element=m.atom_element
                # logger.debug(f'setting mass of {self.name} ({self.type}) to {self.mass} and element to {self.element}')
            else:
                logger.debug(f'no mass record for atom type {self.type} (raw {self.type})')
                exit(-1)
        else:
            logger.error(f'Cannot set mass of {self.name} ({self.type}) from {type(masses)}')
            exit(-1)
    
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

    def set_masses(self,masses):
        for a in self:
            # logger.debug(f'setting mass for {type(a)} \'{str(a)}\'')
            a.set_mass(masses)

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

class CharmmTopDelete:
    def __init__(self,delete_string):
        """
        DELETE ATOM HN ! comment
        """
        self.raw_string=delete_string
        self.error_code=0
        delete_string,self.comment=linesplit(delete_string)
        toks=[x.strip() for x in delete_string.split()]
        if len(toks)<3 or not toks[0].upper().startswith('DELETE'):
            self.error_code=-1
        self.atomname=toks[2].upper()

    def __str__(self):
        return f'DELETE ATOM {self.atomname} ! {self.comment}'

class CharmmTopResi:
    def __init__(self,blockstring,key='RESI',metadata={}):
        self.key=key
        self.blockstring=blockstring
        self.error_code=0
        lines=[x.strip() for x in blockstring.split('\n')]
        # logger.debug(f'processing {key} block with {len(lines)} lines')
        # if len(lines)<2:
        #     logger.debug(f'{lines}')
        #     logger.debug(f'{metadata}')
        titlecard=lines[0]
        assert titlecard.upper().startswith(key),f'bad title card in {key}: [{titlecard}]'
        titledata,titlecomment=linesplit(titlecard)
        tctokens=titledata.split()
        self.resname=tctokens[1]
        self.charge=float(tctokens[2])
        self.synonym=titlecomment.strip()
        self.metadata=metadata

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
        atomgroupcards=list(compress(datacards,isatomgroup))
        # logger.debug(f'{self.resname}: {len(atomgroupcards)} atom group cards')
        # logger.debug(f'{atomgroupcards[0:5]}')
        self.atoms=CharmmTopAtomList([])
        self.atoms_in_group={}
        g=0
        for card in atomgroupcards:
            # logger.debug(f'{card}')
            if card.startswith('ATOM'):
                a=CharmmTopAtom(card)
                self.atoms.append(a)
                if not g in self.atoms_in_group:
                    self.atoms_in_group[g]=CharmmTopAtomList([])
                self.atoms_in_group[g].append(a)
            elif card.startswith('GROU'):
                g+=1
        for a in self.atoms:
            if hasattr(a,'inpatch'):
                if a.inpatch:
                    self.ispatch=True
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
        isDelete=[d.startswith('DELETE') for d in datacards]
        deletecards=compress(datacards,isDelete)
        self.Delete=[]
        for card in deletecards:
            D=CharmmTopDelete(card)
            self.Delete.append(D)

    def set_masses(self,masses):
        self.atoms.set_masses(masses)

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
        logger.debug(f'metadata {m}')
        if m['streamID']=='lipid' and m['substreamID']=='cholesterol':
            self.sterol_annotate()
        elif m['streamID']=='lipid' and m['substreamID'] == 'detergent':
            self.detergent_annotate()
        elif m['streamID']=='lipid' and m['substreamID'] == 'model':
            self.model_annotate()
        else:
            self.generic_lipid_annotate()
    
    def model_annotate(self):
        self.annotation={}
        self.annotation['heads']=[]
        self.annotation['tails']=[]
        self.annotation['shortest_paths']={}

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
        WG=G.copy()
        self.annotation['shortest_paths']={
            head:{
                tail:len(nx.shortest_path(WG,source=head,target=tail)) for tail in self.annotation['tails']
                } for head in self.annotation['heads']
            }

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
