"""

.. module:: residue
   :synopsis: Manages residues
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from .mods import *
# from .config import ConfigGetParam
from .mods import AncestorAwareModList,AncestorAwareMod

class Ter(AncestorAwareMod):
    req_attr=['serial','resname','chainID','resseqnum','insertion']
    yaml_header='Terminals'
    def __init__(self,input_obj):
        if type(input_obj)==dict:
            super().__init__(input_obj)
        elif type(input_obj)==PDBRecord:
            pdbrecord=input_obj
            input_dict={
                'serial':pdbrecord.serial,
                'resname':pdbrecord.residue.resName,
                'chainID':pdbrecord.residue.chainID,
                'resseqnum':pdbrecord.residue.seqNum,
                'insertion':pdbrecord.residue.iCode
            }
            super().__init__(input_dict)
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')
    
class TerList(AncestorAwareModList):
    pass

class Atom(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['serial','name','altloc','resname','chainID','resseqnum','insertion','x','y','z','occ','beta','elem','charge']
    opt_attr=AncestorAwareMod.opt_attr+['segname','empty','link']
    yaml_header='Atoms'
    PDB_keyword='ATOM'

    def __init__(self,input_obj):
        if type(input_obj)==dict:
            super().__init__(input_obj)
        elif type(input_obj)==PDBRecord:
            pdbrecord=input_obj
            input_dict={
                'serial':pdbrecord.serial,
                'name':pdbrecord.name,
                'altloc':pdbrecord.altLoc,
                'resname':pdbrecord.residue.resName,
                'chainID':pdbrecord.residue.chainID,
                'resseqnum':pdbrecord.residue.seqNum,
                'insertion':pdbrecord.residue.iCode,
                'x':pdbrecord.x,
                'y':pdbrecord.y,
                'z':pdbrecord.z,
                'occ':pdbrecord.occupancy,
                'beta':pdbrecord.tempFactor,
                'elem':pdbrecord.element,
                'charge':pdbrecord.charge
            }
            input_dict['segname']=input_dict['chainID']
            input_dict['link']='None'
            input_dict['empty']=False
            super().__init__(input_dict)
        elif type(input_obj)==CIFdict:
            cifdict=input_obj
            al=cifdict['label_alt_id']
            ic=cifdict['pdbx_pdb_ins_code']
            c=cifdict['pdbx_formal_charge']
            input_dict={
                'recordname':'ATOM',
                'serial':int(cifdict['id']),
                'name':cifdict['auth_atom_id'],
                'altloc':' ' if al=='.' else al,
                'resname':cifdict['auth_comp_id'],
                'chainID':cifdict['auth_asym_id'],
                'resseqnum':int(cifdict['auth_seq_id']),
                'insertion':' ' if ic=='?' else ic,
                'x':float(cifdict['cartn_x']),
                'y':float(cifdict['cartn_y']),
                'z':float(cifdict['cartn_z']),
                'occ':float(cifdict['occupancy']),
                'beta':float(cifdict['b_iso_or_equiv']),
                'elem':cifdict['type_symbol'],
                'charge':0.0 if c=='?' else float(c)
            }
            input_dict['segname']=input_dict['chainID']
            input_dict['link']='None'
            input_dict['empty']=False
            super().__init__(input_dict)
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')
    
    def pdb_line(self):
        pdbline='{:<6s}'.format(self.recordname)+\
                '{:5d}'.format(self.serial)+' '+\
                '{:<4s}'.format(' '+self.name if len(self.name)<4 else self.name)+\
                '{:1s}'.format(self.altloc)+\
                '{:<4s}'.format(self.resname)+\
                '{:1s}'.format(self.chainID)+\
                '{:4d}'.format(self.resseqnum)+\
                '{:1s}'.format(self.insertion)+'   '+\
                '{:8.3f}'.format(self.x)+\
                '{:8.3f}'.format(self.y)+\
                '{:8.3f}'.format(self.z)+\
                '{:6.2f}'.format(self.occ)+\
                '{:6.2f}'.format(self.beta)+\
                10*' '+'{:>2s}'.format(self.elem)+'{:2s}'.format(self.charge)
        return pdbline

class AtomList(AncestorAwareModList):
    def adjustSerials(self,Ters):
        ignored_serials=[x.serial for x in Ters]
        if not ignored_serials:
            return
        orig_atom_serials=[x.serial for x in self]
        assert all([not x in orig_atom_serials for x in ignored_serials])
        offending_serial_idx=[]
        cidx=0
        for i in range(len(orig_atom_serials)-1):
            if orig_atom_serials[i]<ignored_serials[cidx]<orig_atom_serials[i+1]:
                offending_serial_idx.append(i+1)
                cidx+=1
                if cidx==len(ignored_serials):
                    break
        for s in offending_serial_idx:
            for i in range(s,len(orig_atom_serials)):
                orig_atom_serials[i]-=1
        for a,s in zip(self,orig_atom_serials):
            if not '_ORIGINAL_' in a.__dict__:
                a._ORIGINAL_={}
            a._ORIGINAL_['serial']=a.serial
            a.serial=s

class Hetatm(Atom):
    PDB_keyword='HETATM'
    yaml_header='Hetatoms'

class Residue(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['resseqnum','insertion','name','chainID','segtype']
    opt_attr=AncestorAwareMod.opt_attr+['atoms','up','down','uplink','downlink','resseqnumi']
    ignore_attr=AncestorAwareMod.ignore_attr+['atoms','up','down','uplink','downlink','resseqnumi']
    _counter=0
    def __init__(self,input_obj):
        if type(input_obj)==dict:
            super().__init__(input_obj)
        elif type(input_obj) in [Atom,Hetatm]:
            a=input_obj
            input_dict={
                'resseqnum':a.resseqnum,
                'insertion':a.insertion,
                'name':a.resname,
                'chainID':a.chainID,
                'atoms':AtomList([a]),
                'segtype':'UNSET',
                'resseqnumi':f'{a.resseqnum}{a.insertion}'
            }
            super().__init__(input_dict)
        elif type(input_obj)==Missing:
            m=input_obj
            input_dict={
                'resseqnum':m.resseqnum,
                'insertion':m.insertion,
                'name':m.resname,
                'chainID':m.chainID,
                'segtype':'UNSET'
            }
            input_dict['resseqnumi']=f'{m.resseqnum}{m.insertion}'
            input_dict['atoms']=AtomList([])
            super().__init__(input_dict)
        else:
            logger.error(f'Cannot initialize {self.__class__} from object type {type(input_obj)}')
        self.index=Residue._counter
        Residue._counter+=1
        self.down=[]
        self.downlink=[]
        self.up=[]
        self.uplink=[]

    def __str__(self):
        return f'{self.chainID}-{self.name}{self.resseqnum}{self.insertion}'
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
        lin=[]
        for d,dl in zip(self.down,self.downlink):
            res.append(d)
            lin.append(dl)
            # logger.debug(f'{str(self)}->{str(d)}')
            tres,tlin=d.get_down_group()
            res.extend(tres)
            lin.extend(tlin)
        return res,lin
    # def resname_charmify(self):
    #     m=ConfigGetParam('PDB_to_CHARMM_Resnames')
    #     if self.name in m:
    #         return m[self.name]
    #     return self.name

class ResidueList(AncestorAwareModList):
    def __init__(self,input_obj):
        if type(input_obj)==list:
            super().__init__(input_obj)
        elif type(input_obj)==AtomList:
            atoms=input_obj
            R=[]
            for a in atoms:
                for r in R[::-1]:
                    if r.add_atom(a):
                        break
                else:
                    R.append(Residue(a))
            super().__init__(R)
        elif type(input_obj)==MissingList:
            R=[Residue(m) for m in input_obj]
            super().__init__(R)

    def get_residue(self,**fields):
        return self.get(**fields)
    def get_atom(self,atname,**fields):
        S=('atoms',{'name':atname})
        return self.get_attr(S,**fields)
    def unique_chainIDs(self):
        c=[]
        for r in self:
            if not r.chainID in c:
                c.append(r.chainID)
        return c
    def atom_serials(self,as_type=str):
        serlist=[]
        for res in self:
            for a in res.atoms:
                serlist.append(as_type(a.serial))
        return serlist
    def caco_str(self,upstream_reslist,seglabel,molid_varname):
        r0=self[0]
        ur=upstream_reslist[-1]
        return f'coord {seglabel} {r0.resseqnum}{r0.insertion} N [cacoIn_nOut {ur.resseqnum}{ur.insertion} {seglabel} ${molid_varname}]'
    def resrange(self,rngrec):
        subR=self.get(chainID=rngrec.chainID)
        subR.sort()
        r1=rngrec.resseqnum1
        i1=rngrec.insertion1
        r2=rngrec.resseqnum2
        i2=rngrec.insertion2
        R1=subR.get(resseqnum=r1,insertion=i1)
        if R1:
            R2=subR.get(resseqnum=r2,insertion=i2)
            if R2:
                idx1=subR.index(R1)
                idx2=subR.index(R2)
                assert idx2>=idx1
                for j in range(idx1,idx2+1):
                    yield self[j]
        return []
    def do_deletions(self,Deletions):
        delete_us=[]
        for d in Deletions:
            for dr in self.resrange(d):
                delete_us.append(dr)
        for d in delete_us:
            self.remove(d)

    def update_links(self,Links:LinkList,atoms:AtomList):
        for link in Links:
            link.residue1=self.get(chainID=link.chainID1,resseqnum=link.resseqnum1,insertion=link.insertion1)
            assert not hasattr(link.residue1,'len'),f'c {link.chainID1} r {link.resseqnum1} i {link.insertion1} returned {len(link.residue1)} residues'
            link.residue2=self.get(chainID=link.chainID2,resseqnum=link.resseqnum2,insertion=link.insertion2)
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
    def apply_segtypes(self,map):
        self.map_attr('segtype','name',map)
