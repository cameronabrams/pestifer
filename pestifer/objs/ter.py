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
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord

from pydantic import Field
from typing import ClassVar
from ..core.baseobj_new import BaseObj, BaseObjList

class Ter(BaseObj):
    """
    A class for handing TER records in PDB files
    """

    _required_fields = {'serial','resname','chainID','resseqnum','insertion'}
    """
    Required attributes for a Ter object.
    These attributes must be provided when creating a Ter object.
    
    - ``serial``: The serial number of the TER record.
    - ``resname``: The residue name of the TER record.
    - ``chainID``: The chain ID of the TER record.
    - ``resseqnum``: The residue sequence number of the TER record.
    - ``insertion``: The insertion code for the TER record.
    """

    serial: int = Field(..., description="Serial number of the TER record")
    resname: str = Field(..., description="Residue name of the TER record")
    chainID: str = Field(..., description="Chain ID of the TER record")
    resseqnum: int = Field(..., description="Residue sequence number of the TER record")
    insertion: str = Field(..., description="Insertion code for the TER record")

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
    
    def describe(self):
        return f"Ter(serial={self.serial}, chainID={self.chainID}, resname={self.resname}, resseqnum={self.resseqnum}, insertion={self.insertion})"
    
    class Adapter:
        def __init__(self, serial: int, resname: str, chainID: str, resseqnum: int, insertion: str):
            self.serial = serial
            self.resname = resname
            self.chainID = chainID
            self.resseqnum = resseqnum
            self.insertion = insertion

        @classmethod
        def from_pdbrecord(cls, pdbrecord: PDBRecord):
            if not pdbrecord.serial:
                # this TER is empty -- probably just a TER on a line by itself in the PDB file
                return cls(serial=-1, resname='', chainID='', resseqnum=-1, insertion='')
            return cls(
                serial=pdbrecord.serial,
                resname=pdbrecord.residue.resName,
                chainID=pdbrecord.residue.chainID,
                resseqnum=pdbrecord.residue.seqNum,
                insertion=pdbrecord.residue.iCode
            )

    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        """
        Create a Ter object from an Adapter instance.
        
        Parameters
        ----------
        adapter : Adapter
            An instance of the Adapter class containing the necessary attributes.
        
        Returns
        -------
        Ter
            A new Ter object created from the attributes of the adapter.
        """
        input_dict = {
            'serial': adapter.serial,
            'resname': adapter.resname,
            'chainID': adapter.chainID,
            'resseqnum': adapter.resseqnum,
            'insertion': adapter.insertion
        }
        return cls(**input_dict)
    
    @singledispatchmethod
    @classmethod
    def new(cls, raw: str) -> "Ter":
        """
        Create a new Ter instance from a shortcode string.
        
        Parameters
        ----------
        raw : str
            The shortcode string in the format serial:resname:chainID:resseqnum:insertion.
        
        Returns
        -------
        Ter
            A new instance of Ter.
        """
        raise TypeError(f"Cannot create Ter from {type(raw)}: {raw}")
    
    @new.register(PDBRecord)
    @classmethod
    def _from_pdbrecord(cls, pdbrecord: PDBRecord) -> "Ter":
        """
        Create a new Ter instance from a PDBRecord.
        
        Parameters
        ----------
        pdbrecord : PDBRecord
            The PDBRecord instance containing the TER information.
        
        Returns
        -------
        Ter
            A new instance of Ter created from the PDBRecord.
        """
        adapter = cls.Adapter.from_pdbrecord(pdbrecord)
        return cls._from_adapter(adapter)
    
class TerList(BaseObjList):
    def describe(self):
        return f'TerList with {len(self)} TER records'
    
    def _validate_item(self, item):
        if not isinstance(item, Ter):
            raise TypeError(f"Expected Ter instance, got {type(item)}")
        return True
