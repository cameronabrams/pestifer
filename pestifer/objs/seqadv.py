# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
``SEQADV`` (PDB) and ``seq_dif`` (mmCIF) records declare differences between the actual sequence of the molecule and the sequence
of its database reference.  This is how authors report things like accidental mutations (conflicts),
engineered mutations, and other differences.  These records are residue-specific.
"""
from __future__ import annotations

import logging

from mmcif.api.PdbxContainers import DataContainer
from pydantic import Field
from pidibble.pdbrecord import PDBRecord, PDBRecordDict
from typing import ClassVar, TYPE_CHECKING

from .resid import ResID
from ..core.baseobj import BaseObj, BaseObjList
from ..util.cifutil import CIFdict

if TYPE_CHECKING:
    from ..molecule.residue import Residue, ResidueList

logger = logging.getLogger(__name__)

class Seqadv(BaseObj):
    """
    A class for handling SEQADV/seq_dif records in input structure files 
    """
    _required_fields = {'idCode', 'resname', 'chainID', 'resid', 'typekey'}
    """
    Required attributes for a Seqadv object.
    These attributes must be provided when creating a Seqadv object.
    
    - ``idCode``: The PDB ID code of the structure.
    - ``resname``: The residue name of the sequence difference.
    - ``chainID``: The chain ID of the segment where the sequence difference occurs.
    - ``resid``: The residue ID where the difference occurs, represented by a `ResID` object.
    - ``typekey``: A key indicating the type of sequence difference (e.g., ``conflict``, ``cloning``, ``expression tag``, ``engineered mutation``, ``variant``, ``insertion``, ``deletion``, ``microheterogeneity``, ``chromophore``, ``user``, ``_other_``).
    """

    _optional_fields = {'database', 'dbAccession', 'dbRes', 'dbSeq', 'pdbx_ordinal', 'pdbx_auth_seq_num', 'residue', 'pdbx_stash'}
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
        'typekey': {'conflict', 'cloning', 'expression tag', 'engineered mutation', 'variant', 'insertion', 'deletion', 'microheterogeneity', 'chromophore', 'user', '_other_'}
    }
    """
    Attribute choices for Seqadv objects.
   
    - ``typekey``: A list of valid sequence difference types, including ``conflict``, ``cloning``, ``expression``, ``engineered``, ``variant``, ``insertion``, ``deletion``, ``microheterogeneity``, ``chromophore``, ``user``, and ``_other_``.
    """

    idCode: str = Field(..., description="PDB ID code of the structure")
    resname: str = Field(..., description="Residue name of the sequence difference")
    chainID: str = Field(..., description="Chain ID of the segment where the sequence difference occurs")
    resid: ResID = Field(..., description="Residue ID where the difference occurs")
    database: str | None = Field(None, description="Database name where the sequence difference is recorded")
    dbAccession: str | None = Field(None, description="Accession number of the sequence in the database")
    dbRes: str | None = Field(None, description="Residue name in the database")
    dbSeq: int | None = Field(None, description="Sequence number in the database")
    typekey: str | None = Field(None, description="Type of sequence difference (e.g., 'conflict', 'cloning', 'expression tag', 'engineered mutation', 'variant', 'insertion', 'deletion', 'microheterogeneity', 'chromophore', 'user', '_other_')")
    pdbx_ordinal: int | None = Field(None, description="Ordinal number of the sequence difference in mmCIF files")
    pdbx_auth_seq_num: int | None = Field(None, description="Author-assigned sequence number in mmCIF files")
    residue: 'Residue' = Field(None, description="Corresponding Residue object, if available")
    pdbx_stash: dict | None = Field(None, description="Stash for temporary storage of attributes")

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

    _CIF_CategoryName: ClassVar[str] = 'struct_ref_seq_dif'
    """
    mmCIF name for Seqadv objects.
    This name is used to identify Seqadv objects in mmCIF files.
    """

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for Seqadv instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        
        Parameters
        ----------
        *args : Any
            The input to adapt, which can be a shortcode string or an Adapter instance.
        
        Returns
        -------
        dict
            A dictionary of parameters for Seqadv instantiation.
        """
        if args and isinstance(args[0], PDBRecord):
            input_dict = Seqadv._from_pdbrecord(args[0])
            return input_dict
        elif args and isinstance(args[0], CIFdict):
            input_dict = Seqadv._from_cifdict(args[0])
            return input_dict
        return super()._adapt(*args, **kwargs)

    @staticmethod
    def _from_pdbrecord(raw: PDBRecord) -> dict:
        """
        Converts a PDBRecord object to a dictionary of parameters for Seqadv.
        """
        if raw.conflict.lower() == 'deletion':
            return dict(
                idCode = raw.idCode,
                chainID = raw.residue.chainID,
                resname = raw.dbRes,
                resid = ResID(raw.dbSeq, raw.dbIcode),
                dbRes = raw.dbRes,
                typekey = 'deletion'
            )
        return dict(
            idCode = raw.idCode,
            chainID = raw.residue.chainID,
            resname = raw.residue.resName,
            resid = ResID(raw.residue.seqNum, raw.residue.iCode),
            dbRes = raw.dbRes,
            typekey = raw.conflict.lower() if raw.conflict.lower() in Seqadv._attr_choices['typekey'] else '_other_',
        )
    
    @staticmethod
    def _from_cifdict(cd: CIFdict) -> dict:
        """
        Converts a CIFdict object to a dictionary of parameters for Seqadv.
        
        Parameters
        ----------
        cd : CIFdict
            The CIFdict instance containing the attributes of the Seqadv object.
        
        Returns
        -------
        dict
            A dictionary of parameters for Seqadv instantiation.
        """
        return {
            'idCode': cd['pdbx_pdb_id_code'],
            'resname': cd['mon_id'],
            'chainID': cd['pdbx_pdb_strand_id'],
            'resid': ResID(cd['pdbx_auth_seq_num'], cd.get('pdbx_auth_seq_ins_code', None)),
            'database': cd.get('pdbx_seq_db_name', None),
            'dbAccession': cd.get('pdbx_seq_db_accession_code', None),
            'dbRes': cd.get('db_mon_id', None),
            'dbSeq': None if cd['pdbx_seq_db_seq_num'] == '' else int(cd['pdbx_seq_db_seq_num']),
            'typekey': cd.get('details', None),
            'pdbx_ordinal': cd.get('pdbx_ordinal', None),
            'pdbx_auth_seq_num': int(cd['pdbx_auth_seq_num'])
        }

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
        text = text.lower()
        keyword_list = self.__class__.attr_choices['typekey']
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
        return f'SEQADV {self.idCode:3s} {self.resname:>3s} {self.chainID:1s} {self.resid.pdbresid:>5s} {self.database:>4s} {self.dbAccession:9s} {self.dbRes:3s} {self.dbSeq:>5d} {self.typekey:21s}          '
    
    def assign_residue(self, Residues: 'ResidueList'):
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
        assert self.residue == None
        if self.typekey != 'deletion':
            if hasattr(self, 'pdbx_auth_seq_num') and self.pdbx_auth_seq_num is not None:
                logger.debug(f'Assigning residue for seqadv {self.chainID}:{self.resid}:{self.typekey} from mmCIF record')
                # this was initialized from a mmCIF record
                # these records are weird in that the 'strand' id
                # is the author-assigned chain id but the sequence
                # number is the mmCIF-assigned sequence number
                self.assign_obj_to_attr('residue', Residues,
                                        auth_asym_id='chainID',
                                        auth_seq_id='pdbx_auth_seq_num')
                if self.residue is not None:
                    if self.pdbx_stash is None:
                        self.pdbx_stash = {}
                    self.pdbx_stash['chainID'] = self.chainID
                    self.pdbx_stash['resid'] = self.resid
                    self.chainID = self.residue.chainID
                    self.resid = self.residue.resid
                    logger.debug(f'...seqadv {self.typekey} auth {self.chainID}:{self.resid} (orig. {self.pdbx_stash["chainID"]}:{self.pdbx_stash["resid"].resid}) assigned to residue {self.residue.resid} (pdb {self.residue.auth_asym_id}:{self.residue.auth_seq_id})')
                else:
                    logger.debug(f'...seqadv {self.typekey} auth {self.chainID}:{self.resid} cannot be resolved from current set of residues')
            else:  # normal PDB record
                # logger.debug(f'Looking to assign residue for seqadv {self.chainID}:{self.resid}')
                self.assign_obj_to_attr('residue', Residues, chainID='chainID', resid='resid')
                if self.residue is None:  # failed to find
                    logger.debug(f'...seqadv {self.typekey} auth {self.chainID}:{self.resid} cannot be resolved from current set of residues')

    def update_from_residue(self):
        """
        Updates the chainID attribute of the seqadv from its residue attribute
        This method sets the chainID attribute of the seqadv to the chainID of its residue
        attribute, if the residue is not None. If the residue is None, it does nothing.
        """
        assert self.residue != None
        self.chainID = self.residue.chainID
        # else:
        #     logger.debug(f'...seqadv {self.typekey} auth {self.pdbx_pdb_strand_id}:{self.pdbx_auth_seq_num} cannot be resolved from current set of residues')
        # we'll assume that if this residue is not found, then this seqadv is never used anyway


class SeqadvList(BaseObjList[Seqadv]):
    """
    A class for handling lists of Seqadvs.
    """

    def describe(self):
        return f'<SeqadvList with {len(self)} seqadvs>'

    @classmethod
    def from_pdb(cls, parsed: PDBRecordDict) -> SeqadvList:
        """
        Create a SeqadvList from a PDBRecordDict.
        """
        if Seqadv._PDB_keyword not in parsed:
            return cls([])
        return cls([Seqadv(x) for x in parsed[Seqadv._PDB_keyword]])

    @classmethod
    def from_cif(cls, parsed: DataContainer) -> SeqadvList:
        """
        Create a SeqadvList from a DataContainer (mmCIF format).
        
        Parameters
        ----------
        parsed : DataContainer
            The parsed mmCIF data container.
        
        Returns
        -------
        SeqadvList
            A new SeqadvList instance containing Seqadv objects created from the mmCIF data.
        """
        obj = parsed.getObj(Seqadv._CIF_CategoryName)
        if obj is not None:
            return cls([Seqadv(CIFdict(obj, i)) for i in range(len(obj))])
        return cls([])

    def assign_residues(self, Residues: 'ResidueList') -> SeqadvList:
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
        delete_us = []
        for s in self.data:
            s.assign_residue(Residues)
        delete_us = self.__class__([s for s in self.data if s.residue is None])
        for s in delete_us:
            self.remove(s)
        return delete_us
    
    def update_from_residues(self):
        """
        Updates the chainID attribute of each Seqadv in the list from its residue.
        This method iterates over each Seqadv in the list and calls its ``update_from_residue`` method.
        It updates the chainID attribute of each Seqadv based on its residue attribute.
        """
        for s in self.data:
            s.update_from_residue()
