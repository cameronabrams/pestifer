# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
``SEQADV`` (PDB) and ``seq_dif`` (mmCIF) records declare differences between the actual sequence of the molecule and the sequence
of its database reference.  This is how authors report things like accidental mutations (conflicts),
engineered mutations, and other differences.  These records are residue-specific.
"""
from pydantic import Field
from typing import ClassVar, Any
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord

from ..core.baseobj_new import BaseObj, BaseObjList

from ..util.cifutil import CIFdict

class Seqadv(BaseObj):
    """
    A class for handling SEQADV/seq_dif records in input structure files 
    """
    _required_fields = ['idCode', 'resname', 'chainID', 'resseqnum', 'insertion', 'typekey']
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
    
    _optional_fields = ['database', 'dbAccession', 'dbRes', 'dbSeq', 'pdbx_ordinal', 'pdbx_auth_seq_num', 'residue']
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

    _attr_choices = {
        'typekey': ['conflict', 'cloning', 'expression', 'engineered', 'variant', 'insertion', 'deletion', 'microheterogeneity', 'chromophore', 'user', '_other_']
    }
    """
    Attribute choices for Seqadv objects.
   
    - ``typekey``: A list of valid sequence difference types, including ``conflict``, ``cloning``, ``expression``, ``engineered``, ``variant``, ``insertion``, ``deletion``, ``microheterogeneity``, ``chromophore``, ``user``, and ``_other_``.
    """

    idCode: str = Field(..., description="PDB ID code of the structure")
    resname: str = Field(..., description="Residue name of the sequence difference")
    chainID: str = Field(..., description="Chain ID of the segment where the sequence difference occurs")
    resseqnum: int = Field(..., description="Residue sequence number where the difference occurs")
    insertion: str = Field(..., description="Insertion code for the residue")
    database: str = Field(None, description="Database name where the sequence difference is recorded")
    dbAccession: str = Field(None, description="Accession number of the sequence in the database")
    dbRes: str = Field(None, description="Residue name in the database")
    dbSeq: int = Field(None, description="Sequence number in the database")
    typekey: str = Field(..., description="Type of sequence difference (e.g., 'conflict', 'cloning', 'expression', 'engineered', 'variant', 'insertion', 'deletion', 'microheterogeneity', 'chromophore', 'user', '_other_')")
    pdbx_ordinal: int = Field(None, description="Ordinal number of the sequence difference in mmCIF files")
    pdbx_auth_seq_num: int = Field(None, description="Author-assigned sequence number in mmCIF files")
    residue: Any = Field(None, description="Corresponding Residue object, if available")

    _yaml_header: ClassVar[str] = 'seqadvs'
    """
    YAML header for Seqadv objects.
    This header is used to identify Seqadv objects in YAML files.
    """

    _objcat: ClassVar[str] = 'seq'
    """
    Object category for Seqadv objects.
    This categorization is used to group Seqadv objects in the object manager.
    """

    _PDB_keyword: ClassVar[str] = 'SEQADV'
    """
    PDB keyword for Seqadv objects.
    This keyword is used to identify Seqadv objects in PDB files.
    """

    _mmCIF_name: ClassVar[str] = 'struct_ref_seq_dif'
    """
    mmCIF name for Seqadv objects.
    This name is used to identify Seqadv objects in mmCIF files.
    """

    def describe(self):
        """
        Describe the Seqadv object.
        
        Returns
        -------
        str
            A string description of the Seqadv object, including idCode, resname, chainID, resseqnum, insertion, and typekey.
        """
        return f"Seqadv(idCode={self.idCode}, resname={self.resname}, chainID={self.chainID}, resseqnum={self.resseqnum}, insertion={self.insertion}, typekey={self.typekey})"
    
    class Adapter:
        def __init__(self, idCode: str, resname: str, chainID: str, resseqnum: int, insertion: str, database: str = None, dbAccession: str = None, dbRes: str = None, dbSeq: int = None, typekey: str = None):
            self.idCode = idCode
            self.resname = resname
            self.chainID = chainID
            self.resseqnum = resseqnum
            self.insertion = insertion
            self.database = database
            self.dbAccession = dbAccession
            self.dbRes = dbRes
            self.dbSeq = dbSeq
            self.typekey = typekey

        @classmethod
        def from_pdbrecord(self, PDBRecord):
            return self(
                idCode=PDBRecord.idCode,
                resname=PDBRecord.residue.resName,
                chainID=PDBRecord.residue.chainID,
                resseqnum=PDBRecord.residue.seqNum,
                insertion=PDBRecord.residue.iCode,
                database=PDBRecord.database,
                dbAccession=PDBRecord.dbAccession,
                dbRes=PDBRecord.dbRes,
                dbSeq=PDBRecord.dbSeq,
                typekey=self.seqadv_details_keyword(PDBRecord.conflict)
            )

        @classmethod
        def from_cifdict(cls, cd: CIFdict):
            input_dict = {
                'idCode': cd['pdbx_pdb_id_code'],
                'resname': cd['mon_id'],
                'chainID': cd['pdbx_pdb_strand_id'],
                'resseqnum': int(cd['seq_num']),
                'insertion': cd['pdbx_pdb_ins_code'],
                'database': cd['pdbx_seq_db_name'],
                'dbAccession': cd['pdbx_seq_db_accession_code'],
                'dbRes': cd['db_mon_id'],
                'dbSeq': int(cd['pdbx_seq_db_seq_num']) if cd['pdbx_seq_db_seq_num'].isdigit() else None,
                'typekey': cls.seqadv_details_keyword(cd['details']),
                'pdbx_auth_seq_num': int(cd['pdbx_auth_seq_num']),
                'pdbx_ordinal': cd.get('pdbx_ordinal', None),
            }
            return cls(**input_dict)

    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        """
        Create a Seqadv object from an Adapter instance, registered by BaseObj.from_input.
        
        Parameters
        ----------
        adapter : Adapter
            The Adapter instance containing the attributes of the Seqadv object.

        Returns
        -------
        Seqadv
            A new Seqadv instance created from the Adapter.
        """
        return cls(**(adapter.__dict__))
    
    @singledispatchmethod
    @classmethod
    def new(cls, raw: Any) -> "Seqadv":
        """
        Create a new Seqadv instance from a shortcode string or other input.
        
        Parameters
        ----------
        raw : str or Adapter
            The shortcode string in the format idCode:resname:chainID:resseqnum-insertion,typekey or an Adapter instance.
        
        Returns
        -------
        Seqadv
            A new instance of Seqadv.
        """
        pass

    @new.register(PDBRecord)
    @classmethod
    def _from_pdbrecord(cls, pdb_record: PDBRecord) -> "Seqadv":
        """
        Create a new Seqadv instance from a PDBRecord.
        
        Parameters
        ----------
        pdb_record : PDBRecord
            The PDBRecord instance containing the attributes of the Seqadv object.
        
        Returns
        -------
        Seqadv
            A new Seqadv instance created from the PDBRecord.
        """
        adapter = cls.Adapter.from_pdbrecord(pdb_record)
        return cls._from_adapter(adapter)
    
    @new.register(CIFdict)
    @classmethod
    def _from_cifdict(cls, cif_dict: CIFdict) -> "Seqadv":
        """
        Create a new Seqadv instance from a CIFdict.
        
        Parameters
        ----------
        cif_dict : CIFdict
            The CIFdict instance containing the attributes of the Seqadv object.
        
        Returns
        -------
        Seqadv
            A new Seqadv instance created from the CIFdict.
        """
        adapter = cls.Adapter.from_cifdict(cif_dict)
        return cls._from_adapter(adapter)

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

class SeqadvList(BaseObjList[Seqadv]):
    """
    A class for handling lists of Seqadvs    
    """

    def describe(self):
        return f'SeqadvList with {len(self)} seqadvs'
    
    def _validate_item(self, item: Seqadv) -> None:
        """
        Validate that the item is an instance of Seqadv.
        
        Parameters
        ----------
        item : Seqadv
            The item to validate.
        
        Raises
        ------
        TypeError
            If the item is not an instance of Seqadv.
        """
        if not isinstance(item, Seqadv):
            raise TypeError(f"Item must be an instance of Seqadv, got {type(item)}")

    def assign_residues(self, Residues):
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
