"""

.. module:: residue
   :synopsis: Manages residues
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from .mods import *
from .config import ConfigGetParam
from .mods import AncestorAwareModList,AncestorAwareMod

class Residue(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['resseqnum','insertion','name','chainID','segtype']
    opt_attr=AncestorAwareMod.opt_attr+['atoms','up','down','uplink','downlink','resseqnumi']
    _counter=0
    def __init__(self,input_dict):
        super().__init__(input_dict)
        self.index=Residue._counter
        Residue._counter+=1
        self.down=[]
        self.downlink=[]
        self.up=[]
        self.uplink=[]

    @classmethod
    def from_atom(cls,a:Atom):
        # print(a.resname,ConfigGetParam('Segtypes_by_Resnames').get(a.resname,'OTHER'))
        input_dict={
            'resseqnum':a.resseqnum,
            'insertion':a.insertion,
            'name':a.resname,
            'chainID':a.chainID,
            'atoms':AtomList([a]),
            'segtype':ConfigGetParam('Segtypes_by_Resnames').get(a.resname,'OTHER'),
            'resseqnumi':f'{a.resseqnum}{a.insertion}'
        }
        inst=cls(input_dict)
        return inst
    @classmethod
    def from_missing(cls,m:Missing):
        input_dict={
            'resseqnum':m.resseqnum,
            'insertion':m.insertion,
            'name':m.resname,
            'chainID':m.chainID,
            'segtype':ConfigGetParam('Segtypes_by_Resnames').get(m.resname,'OTHER')
        }
        input_dict['resseqnumi']=f'{m.resseqnum}{m.insertion}'
        input_dict['atoms']=AtomList([])
        inst=cls(input_dict)
        return inst
    def __str__(self):
        return f'{self.chainID}{self.name}{self.resseqnum}{self.insertion}'
    def __lt__(self,other):
        if self.resseqnum<other.resseqnum:
            return True
        elif self.resseqnum==other.resseqnum:
            if self.insertion==None and other.insertion==None:
                return False
            elif (self.insertion=='' or self.insertion==' ' or self.insertion==None) and other.insertion.isalpha():
                return True
            elif self.insertion.isalpha() and other.insertion.isalpha():
                return ord(self.insertion)<ord(other.insertion)
            else:
                return False
    def add_atom(self,a:Atom):
        if self.resseqnum==a.resseqnum and self.name==a.resname and self.chainID==a.chainID and self.insertion==a.insertion:
            self.atoms.append(a)
            return True
        return False
    def set_chainID(self,chainID):
        self.chainID=chainID
        for a in self.atoms:
            a.chainID=chainID
    def linkTo(self,other,link):
        self.down.append(other)
        self.downlink.append(link)
        other.up.append(self)
        other.uplink.append(link)
    def unlink(self,other,link):
        self.down.remove(other)
        self.downlink.remove(link)
        other.up.remove(self)
        other.uplink.remove(link)
    def ri(self):
        ins0='' if self.insertion==' ' else self.insertion
        return f'{self.resseqnum}{ins0}'
    def get_down_group(self):
        res=[]
        for d in self.down:
            res.append(d)
            res.extend(d.get_down_group())
        return res
    def resname_charmify(self):
        m=ConfigGetParam('PDB_to_CHARMM_Resnames')
        if self.name in m:
            return m[self.name]
        return self.name

class ResidueList(AncestorAwareModList):
    @classmethod
    def from_atoms(cls,atoms:AtomList):
        if not atoms:
            return cls([])
        R=[]
        for a in atoms:
            for r in R[::-1]:
                if r.add_atom(a):
                    break
            else:
                R.append(Residue.from_atom(a))
        return cls(R)
    @classmethod
    def from_missing(cls,missing:MissingList):
        return cls([Residue.from_missing(m) for m in missing])
    def get_residue(self,**fields):
        return self.get(**fields)
    def get_atom(self,atname,**fields):
        S=('atoms',{'name':atname})
        return self.get_attr(S,**fields)
    def atom_serials(self,as_type=str):
        serlist=[]
        for res in self:
            for a in res.atoms:
                serlist.append(as_type(a.serial))
        return serlist
    def caco_str(self,upstream_reslist,seglabel):
        r0=self[0]
        ur=upstream_reslist[-1]
        return 'coord {} {}{} N [cacoIn_nOut {}{} {} 0]'.format(seglabel,r0.resseqnum,r0.insertion,ur.resseqnum,ur.insertion,seglabel)

class LinkList(AncestorAwareModList):
    def update_residues_atoms(self,residues:ResidueList,atoms:AtomList):
        for link in self:
            link.residue1=residues.get(chainID=link.chainID1,resseqnum=link.resseqnum1,insertion=link.iCode1)
            assert not hasattr(link.residue1,'len'),f'c {link.chainID1} r {link.resseqnum1} i {link.iCode} returned {len(link.residue1)} residues'
            link.residue2=residues.get(chainID=link.chainID2,resseqnum=link.resseqnum2,insertion=link.iCode2)
            a1=atoms.get(chainID=link.chainID1,resseqnum=link.resseqnum1,name=link.name1)
            if hasattr(a1,'len'):
                a1=a1.get(altloc=link.altloc1)
            link.atom1=a1
            a2=atoms.get(chainID=link.chainID2,resseqnum=link.resseqnum2,name=link.name2)
            if hasattr(a1,'len'):
                a2=a2.get(altloc=link.altloc2)
            link.atom2=a2
            link.residue1.linkTo(link.residue2,link)
        return self
    def write_TcL(self,transform,mods):
        collect_bytes=''
        for l in self:
            collect_bytes+=l.write_TcL(transform,mods)
        return collect_bytes
