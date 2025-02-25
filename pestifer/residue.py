#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the EmptyResidue and Residue classes
"""
import logging
logger=logging.getLogger(__name__)

from argparse import Namespace
from functools import singledispatchmethod

from pidibble.baserecord import BaseRecord
from pidibble.pdbrecord import PDBRecord

from .atom import Atom, AtomList, Hetatm
from .baseobj import AncestorAwareObj, AncestorAwareObjList
from .config import res_123, segtype_of_resname
from .objs.seqadv import Seqadv, SeqadvList
from .objs.deletion import DeletionList
from .objs.substitution import SubstitutionList
from .stringthings import join_ri, split_ri
from .util.cifutil import CIFdict
from .util.coord import positionN

class EmptyResidue(AncestorAwareObj):
    req_attr=AncestorAwareObj.req_attr+['resname','resseqnum','insertion','chainID','resolved','segtype']
    opt_attr=AncestorAwareObj.opt_attr+['model','id','auth_asym_id','auth_comp_id','auth_seq_id']
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

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # i:C:RRR###A
        # i model number, optional, default 1
        # C chain ID, optional, default 'A'
        # RRR one or three-letter residue code
        # ### integer residue sequence number
        # A insertion code, optional, default ''
        tokens=shortcode.split(':')
        has_modelnum=len(tokens)>2
        model_idx=-1
        chain_idx=-1
        model=1
        chainID='A'
        if has_modelnum:
            model_idx=0
            model=int(tokens[model_idx])
        has_chainID=len(tokens)>1
        if has_chainID:
            chain_idx=model_idx+1
            chainID=tokens[chain_idx]
        res_idx=chain_idx+1
        resncode=tokens[res_idx]
        rescode=''
        ri=''
        okrescode=True
        for byte in resncode:
            if byte.isalpha() and okrescode:
                rescode=rescode+byte
            else:
                okrescode=False
                ri=ri+byte
        rescode=rescode.upper()
        if len(rescode)==1:
            rescode=res_123[rescode]
        r,i=split_ri(ri)
        input_dict={
            'model':model,
            'resname':rescode,
            'chainID':chainID,
            'resseqnum':r,
            'insertion':i,
            'resolved':False,
            'segtype':'UNSET',
        }
        super().__init__(input_dict)

    def __str__(self):
        return f'{self.chainID}_{self.resname}{self.resseqnum}{self.insertion}*'

    def pdb_line(self):
        record_name,code=EmptyResidue.PDB_keyword.split('.')
        return '{:6s}{:>4d}   {:1s} {:3s} {:1s} {:>5d}{:1s}'.format(record_name,
        code,self.model,self.resname,self.chainID,self.resseqnum,self.insertion)

class EmptyResidueList(AncestorAwareObjList):
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
            'model':m.model,
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

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        self.__init__(EmptyResidue(shortcode))

    def __str__(self):
        return super().__str__()[0:-1] # strip off the "*"

    def __lt__(self,other):
        if hasattr(other,'resseqnum') and hasattr(other,'insertion'):
            o_resseqnum=other.resseqnum
            o_insertion=other.insertion
        elif type(other)==str:
            o_resseqnum,o_insertion=split_ri(other)
        else:
            raise Exception(f'I do not know how to make something of type {type(other)} a residue')
        if self.resseqnum<o_resseqnum:
            return True
        elif self.resseqnum==o_resseqnum:
            return self.insertion<o_insertion
        return False
    def __gt__(self,other):
        if hasattr(other,'resseqnum') and hasattr(other,'insertion'):
            o_resseqnum=other.resseqnum
            o_insertion=other.insertion
        elif type(other)==str:
            o_resseqnum,o_insertion=split_ri(other)
        else:
            raise Exception(f'I do not know how to make something of type {type(other)} a residue')
        if self.resseqnum>o_resseqnum:
            return True
        elif self.resseqnum==o_resseqnum:
            return self.insertion>o_insertion
        return False
    def __le__(self,other):
        if self<other:
            return True
        return self.same_resid(other)
    def __ge__(self,other):
        if self>other:
            return True
        return self.same_resid(other)
    
    def same_resid(self,other):
        if type(other)==type(self):
            o_resseqnum=other.resseqnum
            o_insertion=other.insertion
        elif type(other)==str:
            o_resseqnum,o_insertion=split_ri(other)
        return self.resseqnum==o_resseqnum and self.insertion==o_insertion
        
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
        assert type(other)==type(self),f'type of other is {type(other)}; expected {type(self)}'
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
    
class ResidueList(AncestorAwareObjList):

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

    def index(self,R):
        for i,r in enumerate(self):
            if r is R:
                return i
        else:
            raise ValueError(f'Residue not found')

    def map_chainIDs_label_to_auth(self):
        self.chainIDmap_label_to_auth={}
        for r in self:
            if hasattr(r,'auth_asym_id'):
                label_Cid=r.chainID
                auth_Cid=r.auth_asym_id
                if not label_Cid in self.chainIDmap_label_to_auth:
                    self.chainIDmap_label_to_auth[label_Cid]=auth_Cid

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
    def atom_resseqnums(self,as_type=str):
        rlist=[]
        for res in self:
            for a in res.atoms:
                rlist.append(as_type(a.resseqnum))
        return rlist
    def caco_str(self,upstream_reslist,seglabel,molid_varname,tmat):
        r0=self[0]
        ur=upstream_reslist[-1]
        rN=positionN(ur,tmat)
        logger.debug(f'caco {rN}')
        return f'coord {seglabel} {r0.resseqnum}{r0.insertion} N {{{rN[0]:.5f} {rN[1]:.5f} {rN[2]:.5f}}}'
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
        # logger.debug(f'residuelist:apply_segtypes {segtype_of_resname}')
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
            if hasattr(r,'label_asym_id'):
                if not r.chainID in result:
                    result[r.chainID]={}
                if not r.resseqnum in result[r.chainID]:
                    result[r.chainID][r.resseqnum]=Namespace(resseqnum=r.label_seq_id,chainID=r.label_asym_id,insertion=r.insertion)
        return result
    
    def apply_insertions(self,insertions):
        for ins in insertions:
            c,r,i=ins.chainID,ins.resseqnum,ins.insertion
            inc_code=ins.integer_increment
            idx=self.iget(chainID=c,resseqnum=r,insertion=i)
            chainID=self[idx].chainID
            logger.debug(f'insertion begins after {r}{i} which is index {idx} in reslist, chain {chainID}')
            # add residues to residue list
            idx+=1
            i='A' if i in [' ',''] else chr(ord(i)+1)
            for olc in ins.sequence:
                if inc_code:
                    r+=1
                    i=''
                shortcode=f'{chainID}:{olc}{r}{i}'
                new_residue=EmptyResidue(shortcode)
                new_residue.atoms=AtomList([])
                new_residue.segtype='protein'
                self.insert(idx,new_residue)
                logger.debug(f'insertion of new empty residue {shortcode}')
                idx+=1
                i='A' if i in [' ',''] else chr(ord(i)+1)
    
    def renumber(self,links):
        """The possibility exists that empty residues added have resseqnums that conflict with existing resseqnums on the same chain if those resseqnums are in a different segtype (e.g., glycan).  This method will privilege protein residues in such conflicts, and it will renumber non-protein residues, updating any resseqnum records in links """
        protein_residues=self.get(segtype='protein')
        if len(protein_residues)==0: return
        min_protein_resseqnum=min([x.resseqnum for x in protein_residues])
        max_protein_resseqnum=max([x.resseqnum for x in protein_residues])
        non_protein_residues=ResidueList([])
        for p in self:
            if p.segtype!='protein':
                non_protein_residues.append(p)
        logger.debug(f'There are {len(protein_residues)} (resseqnums {min_protein_resseqnum} to {max_protein_resseqnum}) protein residues and {len(non_protein_residues)} non-protein residues')
        assert len(self)==(len(protein_residues)+len(non_protein_residues))
        non_protein_residues_in_conflict=ResidueList([])
        for np in non_protein_residues:
            # logger.debug(f'looking for conflict among protein for chain [{np.chainID}] resseqnum [{np.resseqnum}] insertion [{np.insertion}]')
            tst=protein_residues.get(chainID=np.chainID,resseqnum=np.resseqnum,insertion=np.insertion)
            if tst:
                # logger.debug(f'found it!')
                non_protein_residues_in_conflict.append(np)
        for npc in non_protein_residues_in_conflict:
            non_protein_residues.remove(npc)
        logger.debug(f'There are {len(non_protein_residues_in_conflict)} non-protein residues with resseqnums that conflict with protein residues')

        max_unused_resseqnum=max([max([x.resseqnum for x in protein_residues]),0 if len(non_protein_residues)==0 else max([x.resseqnum for x in non_protein_residues])])+1
        newtst=self.get(resseqnum=max_unused_resseqnum)
        assert newtst==[] #None
        mapper_by_chain={}
        for npc in non_protein_residues_in_conflict:
            c=npc.chainID
            if not c in mapper_by_chain:
                mapper_by_chain[c]={}
            old=join_ri(npc.resseqnum,npc.insertion)
            new=max_unused_resseqnum
            mapper_by_chain[c][old]=new
            max_unused_resseqnum+=1
            npc.resseqnum=new
            npc.insertion=''
            for a in npc.atoms:
                a.resseqnum=new
                a.insertion=''
            logger.debug(f'New resid: {c} {old} -> {new}')
        if mapper_by_chain:
            logger.debug(f'Remapping resids in links')
            for l in links:
                if l.chainID1 in mapper_by_chain:
                    old=join_ri(l.resseqnum1,l.insertion1)
                    if old in mapper_by_chain[l.chainID1]:
                        l.resseqnum1=mapper_by_chain[l.chainID1][old]
                        l.insertion1=''
                        logger.debug(f' remapped 1eft: {l.chainID1} {old} -> {l.resseqnum1}')
                if l.chainID2 in mapper_by_chain:
                    old=join_ri(l.resseqnum2,l.insertion2)
                    if old in mapper_by_chain[l.chainID2]:
                        l.resseqnum2=mapper_by_chain[l.chainID2][old]
                        l.insertion2=''
                        logger.debug(f' remapped right: {l.chainID2} {old} -> {l.resseqnum2}')

    def set_chainIDs(self,chainID):
        for r in self:
            r.set_chainID(chainID)

    def remap_chainIDs(self,the_map):
        for r in self:
            if r.chainID in the_map:
                r.set_chainID(the_map[r.chainID])
