# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..util.cifutil import CIFdict

class Seqadv(AncestorAwareObj):
    """A class for handling SEQADV/seq_dif records in input structure files 
    
    SEQADV/seq_dif records declare differences between the actual sequence of the molecule and the sequence
    of its database reference.  This is how authors report things like accidental mutations (conflicts),
    engineered mutations, and other differences.  These records are residue-specific.

    Attributes
    ----------
    req_attr: list
        * idCode: identification code of the seqadv record
        * resname: residue name
        * chainID: chain identifier
        * resseqnum: sequence number of residue in chain
        * insertion: insertion code
        * typekey: reports the type of the seqadv/seq_dif based on the description provided by the author
    opt_attr: list
        * database: name of database
        * dbAccession: accession code of this molecule's sequence in that database
        * dbRes: name of residue in the database
        * dbSeq: sequence position of residue in the database
        * pdbx_pdb_strand_id: (mmCIF format) chain identifier according to author
        * pdbx_auth_seq_num: (mmCIF format) residue sequence position according to author
        * pdbx_ordinal: some integer, i don't know
        * residue: holder for actual residue this seqadv refers to

        _from_cifdict   _struct_ref_seq_dif.align_id 
                        _struct_ref_seq_dif.pdbx_pdb_id_code 
                        _struct_ref_seq_dif.mon_id 
                        _struct_ref_seq_dif.pdbx_pdb_strand_id 
                        _struct_ref_seq_dif.seq_num 
                        _struct_ref_seq_dif.pdbx_pdb_ins_code 
                        _struct_ref_seq_dif.pdbx_seq_db_name 
                        _struct_ref_seq_dif.pdbx_seq_db_accession_code 
                        _struct_ref_seq_dif.db_mon_id 
                        _struct_ref_seq_dif.pdbx_seq_db_seq_num 
                        _struct_ref_seq_dif.details 
                        _struct_ref_seq_dif.pdbx_auth_seq_num 
                        _struct_ref_seq_dif.pdbx_ordinal 

    Methods
    -------
    seqadv_details_keyword(text)
        returns the typekey for the given text
    
    pdb_line()
        write the seqadv as it would appear in a PDB file
    
    assign_residue(Residues)
        given the list of residues 'Residues', finds
        the residue that matches the specs of the caller
        assigns it to the 'residue ' element
    
    update_residue()
        assigns the chainID attribute of the caller to
        the value of the chainID of the caller's 
        residue attribute
    
    """
    req_attr=AncestorAwareObj.req_attr+['idCode','resname','chainID','resseqnum','insertion','typekey']
    opt_attr=AncestorAwareObj.opt_attr+['database','dbAccession','dbRes','dbSeq','pdbx_ordinal','pdbx_auth_seq_num','residue']
    attr_choices=AncestorAwareObj.attr_choices.copy()
    attr_choices.update({'typekey':['conflict','cloning','expression','typekey','engineered','variant','insertion','deletion','microheterogeneity','chromophore','user','_other_']})
    yaml_header='seqadvs'
    objcat='seq'
    PDB_keyword='SEQADV'
    mmCIF_name='struct_ref_seq_dif'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(PDBRecord)
    def _from_pdbrecord(self,pdbrecord):
        input_dict={
            'idCode':pdbrecord.idCode,
            'resname':pdbrecord.residue.resName,
            'chainID':pdbrecord.residue.chainID,
            'resseqnum':pdbrecord.residue.seqNum,
            'insertion':pdbrecord.residue.iCode,
            'database':pdbrecord.database,
            'dbAccession':pdbrecord.dbAccession,
            'dbRes':pdbrecord.dbRes,
            'dbSeq':pdbrecord.dbSeq,
            'typekey':self.seqadv_details_keyword(pdbrecord.conflict),
            'residue':None
        }
        super().__init__(input_dict)

    @__init__.register(CIFdict)
    def _from_cifdict(self,cd):

        input_dict={
            'idCode':cd['pdbx_pdb_id_code'],
            'resname':cd['mon_id'],
            'chainID':cd['pdbx_pdb_strand_id'],
            'resseqnum':cd['seq_num'],
            'insertion':cd['pdbx_pdb_ins_code'],
            'database':cd['pdbx_seq_db_name'],
            'dbAccession':cd['pdbx_seq_db_accession_code'],
            'dbRes':cd['db_mon_id'],
            'dbSeq':cd['pdbx_seq_db_seq_num'],
            'typekey':self.seqadv_details_keyword(cd['details']),
            'pdbx_auth_seq_num':cd['pdbx_auth_seq_num'],
            'pdbx_ordinal':cd['pdbx_ordinal'],
            'label_asym_id':'UNSET',
            'residue':None
        }
        input_dict['resseqnum']=int(input_dict['resseqnum'])
        input_dict['pdbx_auth_seq_num']=int(input_dict['pdbx_auth_seq_num'])
        if input_dict['dbSeq'].isdigit():
            input_dict['dbSeq']=int(input_dict['dbSeq'])
        super().__init__(input_dict)

    def seqadv_details_keyword(self,text):
        text=text.lower()
        keyword_list=self.__class__.attr_choices['typekey']
        for keyword in keyword_list:
            if keyword in text:
                return keyword
        return '_other_'

    def pdb_line(self):
        return f'SEQADV {self.idCode:3s} {self.resname:>3s} {self.chainID:1s} {self.resseqnum:>4d}{self.insertion:1s} {self.database:>4s} {self.dbAccession:9s} {self.dbRes:3s} {self.dbSeq:>5d} {self.typekey:21s}          '
    
    def assign_residue(self,Residues):
        # logger.debug(f'Searching {len(Residues)} Residues for auth_chain {self.pdbx_pdb_strand_id} auth_seq {self.pdbx_auth_seq_num}')
        assert self.residue==None
        if self.typekey!='deletion':
            if hasattr(self,'pdbx_auth_seq_num'):
                # this was initialized from a mmCIF record
                # these records are weird in that the 'strand' id
                # is the author-assigned chain id but the sequence
                # number is the mmCIF-assigned sequence number
                self.assign_obj_to_attr('residue',Residues,
                                auth_asym_id='chainID',
                                auth_seq_id='pdbx_auth_seq_num',
                                insertion='insertion')
                if self.residue==[]:
                    self.residue=None
                if self.residue!=None:
                    self.chainID=self.residue.chainID
            else: # normal PDB record
                self.assign_obj_to_attr('residue',Residues,chainID='chainID',resseqnum='resseqnum',insertion='insertion')
                if self.residue==[]: # failed to find 
                    self.residue=None

    def update_from_residue(self):
        assert self.residue!=None
        self.chainID=self.residue.chainID
        # else:
        #     logger.debug(f'...seqadv {self.typekey} auth {self.pdbx_pdb_strand_id}:{self.pdbx_auth_seq_num} cannot be resolved from current set of residues')
        # we'll assume that if this residue is not found, then this seqadv is never used anyway

class SeqadvList(AncestorAwareObjList):
    """A class for handling lists of Seqadvs
    
    Methods
    -------
    assign_residues(Residues):
        calls assign_residue for each member element
        and returns a list of residues that were
        not assigned (out of the list Residues)
    
    update_from_residues()
        calls update_from_residue for each member element
    
    """
    def assign_residues(self,Residues):
        delete_us=[]
        for s in self:
            s.assign_residue(Residues)
        delete_us=self.__class__([s for s in self if s.residue==None])
        for s in delete_us:
            self.remove(s)
        return delete_us
    
    def update_from_residues(self):
        for s in self:
            s.update_from_residue()
