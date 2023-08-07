"""

.. module:: residue
   :synopsis: Manages residues
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from .mods import *
from .config import ConfigGetParam

class Residue(CloneableMod):
    req_attr=['resseqnum','insertion','name','chainID','segtype']
    opt_attr=['atoms','up','down','uplink','downlink']
    _counter=0
    def __init__(self,input_dict):
        super().__init__(input_dict)
        self.index=Residue._counter
        Residue._counter+=1
        self.segtype=''
        self.down=[]
        self.downlink=[]
        self.up=[]
        self.uplink=[]

    @classmethod
    def from_atom(cls,a:Atom):
        input_dict={
            'resseqnum':a.resseqnum,
            'insertion':a.insertion,
            'name':a.resname,
            'chainID':a.chainID,
            'atoms':AtomList([a]),
            'segtype':ConfigGetParam('Segtypes_by_Resnames').get(a.resname,'OTHER')
        }
        input_dict['resseqnumi']=f'{a.resseqnum}{a.insertion}'
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
        return self.nam

class ResidueList(CloneableModList):
    def get_residue(self,**fields):
        return self.get(**fields)
    def get_atom(self,atname,**fields):
        S=('atoms',{'name':atname})
        return self.get_attr(S,**fields)
    def atom_serials(self):
        serlist=[]
        for a in self.atoms:
            serlist.append(a.serial)
        return serlist
    