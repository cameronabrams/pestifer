# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
``SEQADV`` (PDB) and ``seq_dif`` (mmCIF) records declare differences between the actual sequence of the molecule and the sequence
of its database reference.  This is how authors report things like accidental mutations (conflicts),
engineered mutations, and other differences.  These records are residue-specific.
"""
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..util.cifutil import CIFdict

class Seqadv(AncestorAwareObj):
    """
    A class for handling SEQADV/seq_dif records in input structure files 
    """

    req_attr=AncestorAwareObj.req_attr+['idCode','resname','chainID','resseqnum','insertion','typekey']
    """
    Required attributes for a Seqadv object.
    These attributes must be provided when creating a Seqadv object.
    
    - ``idCode``: The PDB ID code of the structure.
    - ``resname``: The residue name of the sequence difference.
    - ``chainID``: The chain ID of the segment where the sequence difference occurs.
    - ``resseqnum``: The residue sequence number where the difference occurs.
    - ``insertion``: The insertion code for the residue.
    - ``typekey``: A key indicating the type of sequence difference (e.g., ``conflict``, ``cloning``, ``expression``, ``engineered``, ``variant``, ``insertion``, ``deletion``, ``microheterogeneity``, ``chromophore``, ``user``, ``_other_``).
    """
    
    opt_attr=AncestorAwareObj.opt_attr+['database','dbAccession','dbRes','dbSeq','pdbx_ordinal','pdbx_auth_seq_num','residue']
    """
    Optional attributes for a Seqadv object.
    These attributes may be present but are not required.
    
    - ``database``: The database name where the sequence difference is recorded.
    - ``dbAccession``: The accession number of the sequence in the database.
    - ``dbRes``: The residue name in the database.
    - ``dbSeq``: The sequence number in the database.
    - ``pdbx_ordinal``: The ordinal number of the sequence difference in mmCIF files.
    - ``pdbx_auth_seq_num``: The author-assigned sequence number in mmCIF files.
    - ``residue``: The corresponding Residue object, if available. This attribute is used to link the Seqadv object to a specific Residue object in the structure, allowing for easier access to the residue's properties.
    """
    attr_choices=AncestorAwareObj.attr_choices.copy()
    """
    Attribute choices for Seqadv objects.
   
    - ``typekey``: A list of valid sequence difference types, including ``conflict``, ``cloning``, ``expression``, ``engineered``, ``variant``, ``insertion``, ``deletion``, ``microheterogeneity``, ``chromophore``, ``user``, and ``_other_``.
    """
    attr_choices.update({'typekey':['conflict','cloning','expression','typekey','engineered','variant','insertion','deletion','microheterogeneity','chromophore','user','_other_']})

    yaml_header='seqadvs'
    """
    YAML header for Seqadv objects.
    This header is used to identify Seqadv objects in YAML files.
    """
    
    objcat='seq'
    """
    Object category for Seqadv objects.
    This categorization is used to group Seqadv objects in the object manager.
    """

    PDB_keyword='SEQADV'
    """
    PDB keyword for Seqadv objects.
    This keyword is used to identify Seqadv objects in PDB files.
    """
    
    mmCIF_name='struct_ref_seq_dif'
    """
    mmCIF name for Seqadv objects.
    This name is used to identify Seqadv objects in mmCIF files.
    """
    
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
        """
        Returns the typekey for the given text
        
        Parameters
        ----------
        text : str
            The text to analyze for keywords.
        
        Returns
        -------
        str
            The typekey corresponding to the text, or ``_other_`` if no keyword is found.
        """
        text=text.lower()
        keyword_list=self.__class__.attr_choices['typekey']
        for keyword in keyword_list:
            if keyword in text:
                return keyword
        return '_other_'

    def pdb_line(self):
        """
        Write the seqadv as it would appear in a PDB file
        
        Returns
        -------
        str
            The PDB line for the seqadv record.
        """
        return f'SEQADV {self.idCode:3s} {self.resname:>3s} {self.chainID:1s} {self.resseqnum:>4d}{self.insertion:1s} {self.database:>4s} {self.dbAccession:9s} {self.dbRes:3s} {self.dbSeq:>5d} {self.typekey:21s}          '
    
    def assign_residue(self,Residues):
        """
        Assigns the residue attribute of the seqadv to the corresponding Residue object
        This method searches through the provided list of Residues to find a match based on the
        attributes of the seqadv object. If a match is found, the residue attribute is set
        to the corresponding Residue object. If no match is found, the residue attribute remains None.
        
        Parameters
        ----------
        Residues : list
            A list of Residue objects to search for the corresponding residue.
        """
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
        """
        Updates the chainID attribute of the seqadv from its residue attribute
        This method sets the chainID attribute of the seqadv to the chainID of its residue
        attribute, if the residue is not None. If the residue is None, it does nothing.
        """
        assert self.residue!=None
        self.chainID=self.residue.chainID
        # else:
        #     logger.debug(f'...seqadv {self.typekey} auth {self.pdbx_pdb_strand_id}:{self.pdbx_auth_seq_num} cannot be resolved from current set of residues')
        # we'll assume that if this residue is not found, then this seqadv is never used anyway

class SeqadvList(AncestorAwareObjList):
    """
    A class for handling lists of Seqadvs    
    """
    def assign_residues(self,Residues):
        """
        Assigns residues to each Seqadv in the list from the provided Residues. This method iterates over each Seqadv in the list and calls its ``assign_residue`` method  with the provided Residues. After assigning residues, it creates a new SeqadvList containing
        the Seqadv objects that have no assigned residue (i.e., their ``residue`` attribute is None).
        It then removes these Seqadv objects from the original list.
        
        Parameters
        ----------
        Residues : list
            A list of Residue objects to assign to the Seqadv objects.
        
        Returns
        -------
        SeqadvList
            A new SeqadvList containing the Seqadv objects that could not be assigned a residue.

        """        
        delete_us=[]
        for s in self:
            s.assign_residue(Residues)
        delete_us=self.__class__([s for s in self if s.residue==None])
        for s in delete_us:
            self.remove(s)
        return delete_us
    
    def update_from_residues(self):
        """
        Updates the chainID attribute of each Seqadv in the list from its residue.
        This method iterates over each Seqadv in the list and calls its ``update_from_residue`` method.
        It updates the chainID attribute of each Seqadv based on its residue attribute.
        """
        for s in self:
            s.update_from_residue()
