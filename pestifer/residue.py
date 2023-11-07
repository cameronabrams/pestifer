#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Atoms and residues
"""
from .mods import *
from .config import segtype_of_resname
from pidibble.baserecord import BaseRecord
from functools import singledispatchmethod
from argparse import Namespace

class Atom(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['serial','name','altloc','resname','chainID','resseqnum','insertion','x','y','z','occ','beta','elem','charge']
    opt_attr=AncestorAwareMod.opt_attr+['segname','empty','link','recordname','auth_seq_id','auth_comp_id','auth_asym_id','auth_atom_id']
    yaml_header='atoms'
    PDB_keyword='ATOM'
    mmCIF_name='atom_site'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(BaseRecord)
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
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

    @__init__.register(CIFdict)
    def _from_cifdict(self,cifdict):
        seq_id=cifdict['label_seq_id']
        if seq_id=='':
            seq_id=cifdict['auth_seq_id']
        input_dict={
            'recordname':'ATOM',
            'serial':int(cifdict['id']),
            'name':cifdict['label_atom_id'],
            'altloc':cifdict['label_alt_id'],
            'resname':cifdict['label_comp_id'],
            'chainID':cifdict['label_asym_id'],
            'resseqnum':int(seq_id),
            'insertion':cifdict['pdbx_pdb_ins_code'],
            'x':float(cifdict['cartn_x']),
            'y':float(cifdict['cartn_y']),
            'z':float(cifdict['cartn_z']),
            'occ':float(cifdict['occupancy']),
            'beta':float(cifdict['b_iso_or_equiv']),
            'elem':cifdict['type_symbol'],
            'charge':cifdict['pdbx_formal_charge'],
            'auth_atom_id':cifdict['auth_atom_id'],
            'auth_seq_id':int(cifdict['auth_seq_id']),
            'auth_comp_id':cifdict['auth_comp_id'],
            'auth_asym_id':cifdict['auth_asym_id']
        }
        input_dict['segname']=input_dict['chainID']
        input_dict['link']='None'
        input_dict['empty']=False
        super().__init__(input_dict)
    
    def pdb_line(self):
        pdbline='{:<6s}'.format(self.recordname)+\
                '{:5d}'.format(self.serial)+' '+\
                '{:<4s}'.format(' '+self.resname if len(self.resname)<4 else self.resname)+\
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
    yaml_header='hetatoms'

class EmptyResidue(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['resname','resseqnum','insertion','chainID','resolved','segtype']
    opt_attr=AncestorAwareMod.opt_attr+['model','id','auth_asym_id','auth_comp_id','auth_seq_id']
    yaml_header='missings'
    PDB_keyword='REMARK.465'
    mmCIF_name='pdbx_unobs_or_zero_occ_residues'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(BaseRecord)
    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'model':pdbrecord.modelNum,
            'resname':pdbrecord.resName,
            'chainID':pdbrecord.chainID,
            'resseqnum':pdbrecord.seqNum,
            'insertion':pdbrecord.iCode,
            'resolved':False,
            'segtype':'UNSET'
        }
        super().__init__(input_dict)
    
    @__init__.register(CIFdict)
    def _from_cifdict(self,cd):
        mn=cd['pdb_model_num']
        if type(mn)==str and mn.isdigit:
            nmn=int(mn)
        else:
            nmn=1
        input_dict={
            'model':nmn,
            'resname':cd['label_comp_id'],
            'chainID':cd['label_asym_id'],
            'resseqnum':int(cd['label_seq_id']),
            'insertion':cd['pdb_ins_code'],
            'resolved':False,
            'segtype':'UNSET',
            'auth_asym_id':cd['auth_asym_id'],
            'auth_comp_id':cd['auth_comp_id'],
            'auth_seq_id':int(cd['auth_seq_id']),
        }
        super().__init__(input_dict)

    def pdb_line(self):
        record_name,code=EmptyResidue.PDB_keyword.split('.')
        return '{:6s}{:>4d}   {:1s} {:3s} {:1s} {:>5d}{:1s}'.format(record_name,
        code,self.model,self.resname,self.chainID,self.resseqnum,self.insertion)

class EmptyResidueList(AncestorAwareModList):
    pass

class Residue(EmptyResidue):
    req_attr=EmptyResidue.req_attr+['atoms']
    opt_attr=EmptyResidue.opt_attr+['up','down','uplink','downlink']
    ignore_attr=EmptyResidue.ignore_attr+['atoms','up','down','uplink','downlink']
    _counter=0

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    
    @__init__.register(Atom)
    @__init__.register(Hetatm)
    def _from_atom(self,a):
        input_dict={
            'resseqnum':a.resseqnum,
            'insertion':a.insertion,
            'resname':a.resname,
            'chainID':a.chainID,
            'atoms':AtomList([a]),
            'segtype':'UNSET',
            'resolved':True
        }
        for cif_xtra in ['auth_seq_id','auth_comp_id','auth_asym_id']:
            if hasattr(a,cif_xtra):
                input_dict[cif_xtra]=a.__dict__[cif_xtra]
        super().__init__(input_dict)
        self.index=Residue._counter
        Residue._counter+=1
        self.down=[]
        self.downlink=[]
        self.up=[]
        self.uplink=[]

    @__init__.register(EmptyResidue)
    def _from_emptyresidue(self,m):
        input_dict={
            'resseqnum':m.resseqnum,
            'insertion':m.insertion,
            'resname':m.resname,
            'chainID':m.chainID,
            'atoms':AtomList([]),
            'segtype':'UNSET',
            'resolved':False
        }
        for cif_xtra in ['auth_asym_id','auth_comp_id','auth_seq_id']:
            if hasattr(m,cif_xtra):
                input_dict[cif_xtra]=m.__dict__[cif_xtra]
        super().__init__(input_dict)
        self.index=Residue._counter
        Residue._counter+=1
        self.down=[]
        self.downlink=[]
        self.up=[]
        self.uplink=[]

    def __str__(self):
        rc='' if self.resolved else '*'
        return f'{self.chainID}-{self.resname}{self.resseqnum}{self.insertion}{rc}'
    def __lt__(self,other):
        if self.resseqnum<other.resseqnum:
            return True
        elif self.resseqnum==other.resseqnum:
            return self.insertion<other.insertion
        return False
    def __le__(self,other):
        if self<other:
            return True
        return self.resseqnum==other.resseqnum and self.insertion==other.insertion
        
    def add_atom(self,a:Atom):
        if self.resseqnum==a.resseqnum and self.resname==a.resname and self.chainID==a.chainID and self.insertion==a.insertion:
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
    
class ResidueList(AncestorAwareModList):

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(AtomList)
    def _from_atomlist(self,atoms):
        R=[]
        for a in atoms:
            for r in R[::-1]:
                if r.add_atom(a):
                    break
            else:
                R.append(Residue(a))
        super().__init__(R)
    
    @__init__.register(EmptyResidueList)
    def _from_emptyresiduelist(self,input_list):
        R=[Residue(m) for m in input_list]
        super().__init__(R)

    def map_chainIDs_label_to_auth(self):
        self.chainIDmap_cif_to_pdb={}
        for r in self:
            if hasattr(r,'auth_asym_id'):
                aCid=r.auth_asym_id
                Cid=r.chainID
                if not Cid in self.chainIDmap_cif_to_pdb:
                    self.chainIDmap_cif_to_pdb[Cid]=aCid
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

    def apply_segtypes(self):
        self.map_attr('segtype','resname',segtype_of_resname)
    
    def deletion(self,DL:DeletionList):
        excised=[]
        for d in DL:
            chain=self.get(chainID=d.chainID)
            r1=chain.get(resseqnum=d.resseqnum1,insertion=d.insertion1)
            r2=chain.get(resseqnum=d.resseqnum2,insertion=d.insertion2)
            for r in chain:
                if r1<=r<=r2:
                    excised.append(r)
        for x in excised:
            self.remove(x)
            assert not x in self
        return excised

    def substitutions(self,SL:SubstitutionList):
        delete_us=[]
        newseqadv=SeqadvList([]) # for holding single-residue changes for resolved residues
        for s in SL:
            subseq=s.subseq
            currsubidx=0
            chain=self.get(chainID=s.chainID)
            assert chain!=None,f'Error: no chain {s.chainID}'
            r1=chain.get(resseqnum=s.resseqnum1,insertion=s.insertion1)
            assert r1!=None,f'Error: no resseqnum {s.resseqnum1} insertion [{s.insertion1}]'
            r2=chain.get(resseqnum=s.resseqnum2,insertion=s.insertion2)
            assert r2!=None,f'Error: no resseqnum {s.resseqnum2} insertion [{s.insertion2}]'
            for r in chain:
                if r1<=r<=r2:
                    if currsubidx<len(subseq):
                        resname=res_123[subseq[currsubidx].upper()]
                        if r.resolved: # make a new seqadv for this mutation
                            input_dict={
                                'idCode':'I doubt I ever use this',
                                'typekey':'user',
                                'resname':r.resname,
                                'chainID':r.chainID,
                                'resseqnum':r.resseqnum,
                                'insertion':r.insertion,
                                'dbRes':resname
                            }
                            newseqadv.append(Seqadv(input_dict))
                        else:  # just change the residue name
                            r.name=resname
                        currsubidx+=1
                    else:
                        delete_us.append(r)
            if currsubidx<len(subseq):
                # we have unsubsituted residue(s) left that must be inserted
                pass
        for r in delete_us:
            self.remove(r)
        return newseqadv,delete_us

    def cif_residue_map(self):
        result={}
        for r in self:
            if hasattr(r,'auth_asym_id'):
                if not r.chainID in result:
                    result[r.chainID]={}
                if not r.resseqnum in result[r.chainID]:
                    result[r.chainID][r.resseqnum]=Namespace(resseqnum=r.auth_seq_id,chainID=r.auth_asym_id,insertion=r.insertion)
        return result