# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for handling TER records in PDB files.     
Why is this necessary?  TER records are used to indicate
'breaks' in the list of ATOM records in a PDB file to signify
physical discontinuity of the amino acid chain(s).  The problem
we have to solve is that sometimes a TER record gets a unique 
atom serial number, even though they are not atoms. So in a list 
of atoms the sequence of serial numbers can be discontinuous.
This is not a problem per se, *but* when VMD reads in a PDB file 
it ignores both TER records *and* the explicit atom serial numbers.
Pidibble by default does not, so to make the atom records 
VMD-compliant, we have to adjust atom serial numbers, and to do 
that we have to know where any non-empty TER records are.
"""

import logging

from pidibble.pdbrecord import PDBRecord, PDBRecordDict
from pydantic import Field
from typing import ClassVar

from .resid import ResID

from ..core.baseobj import BaseObj, BaseObjList

logger = logging.getLogger(__name__)

class Ter(BaseObj):
    """
    A class for handing TER records in PDB files
    """

    _optional_fields = {'serial', 'resname', 'chainID', 'resid'}
    """
    Optional attributes for a Ter object.  Ters can be empty.
    These attributes can be provided when creating a Ter object.
    
    - ``serial``: The serial number of the TER record.
    - ``resname``: The residue name of the TER record.
    - ``chainID``: The chain ID of the TER record.
    - ``resid``: The residue ID of the TER record.
    """

    serial: int | None = Field(None, description="Serial number of the TER record")
    resname: str | None = Field(None, description="Residue name of the TER record")
    chainID: str | None = Field(None, description="Chain ID of the TER record")
    resid: ResID | None = Field(None, description="Residue ID of the TER record")

    _yaml_header: ClassVar[str] ='terminals'
    """ 
    YAML header for Ter objects.
    This header is used to identify Ter objects in YAML files.
    """

    _objcat: ClassVar[str] = 'seq'
    """
    Category of the Ter object.
    This categorization is used to group Ter objects in the object manager.
    """

    _PDB_keyword: ClassVar[str] = 'TER'
    """
    The PDB keyword for TER records.
    This keyword is used to identify TER records in PDB files.
    """
    
    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for Ter instantiations.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if args and isinstance(args[0], PDBRecord):
            pdbrecord: PDBRecord = args[0]
            assert pdbrecord.key == 'TER'
            # logger.debug(f'serial "{pdbrecord.serial}" ({type(pdbrecord.serial)})')
            # logger.debug(f'Adapting PDBRecord {str(pdbrecord)} to Ter')
            if not pdbrecord.serial:
                # this TER is empty
                # logger.debug(f'Empty ter')
                input_dict = {
                    'serial': None,
                    'resname': None,
                    'chainID': None,
                    'resid': None
                }
                return input_dict
            else:
                input_dict = {
                    'serial': pdbrecord.serial,
                    'resname': pdbrecord.residue.resName,
                    'chainID': pdbrecord.residue.chainID,
                    'resid': ResID(resseqnum=pdbrecord.residue.seqNum, insertion=pdbrecord.residue.iCode)
                }
            return input_dict
        return super()._adapt(*args, **kwargs)

class TerList(BaseObjList[Ter]):

    # _has_serials: ClassVar[bool] = False

    def describe(self):
        return f'<TerList with {len(self)} TER records, has_serials={self._has_serials}>'

    @classmethod
    def from_pdb(cls, parsed: PDBRecordDict, model_id = None) -> "TerList":
        """
        Create a TerList from parsed PDB data.
        
        Parameters
        ----------
        parsed : PDBRecordDict
            A dictionary containing parsed PDB records.
        model_id : int, optional
            The model ID to filter the TER records by. If None, all models are included.
        
        Returns
        -------
        TerList
            A new TerList instance containing Ter objects created from the PDB data.
        """
        if Ter._PDB_keyword not in parsed or not parsed[Ter._PDB_keyword]:
            return cls([])
        # logger.debug(f'Reading from {[repr(p) for p in parsed[Ter._PDB_keyword]]}')
        L = cls([Ter(p) for p in parsed[Ter._PDB_keyword] if (model_id is None or p.model == model_id)])
        return L

    def has_serials(self):
        return any(x.serial is not None for x in self)