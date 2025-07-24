# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
CleavageSite is a class for handling chain cleavage in molecular structures.
It represents a user-specified modification that cleaves a segment of a chain
between two specified residues, creating a new segment.
"""
import logging
logger=logging.getLogger(__name__)

from pydantic import Field
from typing import ClassVar
from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.stringthings import split_ri, join_ri

class CleavageSite(BaseObj):
    """
    A class for handling chain cleavage.  Note that this mod is not expected to be part of an 
    ObjManager so the yaml_header and objcat attributes are irrelevant.  It is instead handled
    as a run task.
    """
    _required_fields = ['chainID', 'resseqnum1', 'insertion1', 'resseqnum2', 'insertion2']

    chainID:    str = Field(..., description="Chain ID of the segment to be cleaved")
    resseqnum1: int = Field(..., description="N-terminal residue number of the cleavage site")
    insertion1: str = Field(..., description="Insertion code of the N-terminal residue")
    resseqnum2: int = Field(..., description="C-terminal residue number of the cleavage site")
    insertion2: str = Field(..., description="Insertion code of the C-terminal residue")

    _yaml_header: ClassVar[str] = 'cleavages'
    _objcat: ClassVar[str] = 'seq'

    def describe(self):
        return f"CleavageSite(chainID={self.chainID}, resseqnum1={self.resseqnum1}, insertion1={self.insertion1}, resseqnum2={self.resseqnum2}, insertion2={self.insertion2})"

    class Adapter:
        """
        A class to represent the shortcode format for CleavageSite, so that we can register to BaseObj.from_input rather than defining a local from_input.

        The shortcode format is C:R1-R2
        where:
        - C is the chain ID
        - R1 is the residue number and insertion code of the N-terminal partner in the peptide bond
        - R2 is the residue number and insertion code of the C-terminal partner in the peptide bond
        """
        def __init__(self, chainID: str, resseqnum1: int, insertion1: str, resseqnum2: int, insertion2: str):
            self.chainID = chainID
            self.resseqnum1 = resseqnum1
            self.insertion1 = insertion1
            self.resseqnum2 = resseqnum2
            self.insertion2 = insertion2

        @classmethod
        def from_string(cls, raw: str):
            chainID, res_range = raw.split(":")
            resseqnum1, insertion1 = split_ri(res_range.split("-")[0])
            resseqnum2, insertion2 = split_ri(res_range.split("-")[1])
            return cls(chainID, resseqnum1, insertion1, resseqnum2, insertion2)

        def to_string(self) -> str:
            return f"{self.chainID}:{join_ri(self.resseqnum1, self.insertion1)}-{join_ri(self.resseqnum2, self.insertion2)}"


    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_shortcode(cls,shortcode:Adapter):
        input_dict={
            'chainID':shortcode.chainID,
            'resseqnum1':shortcode.resseqnum1,
            'insertion1':shortcode.insertion1,
            'resseqnum2':shortcode.resseqnum2,
            'insertion2':shortcode.insertion2
        }
        return cls(**input_dict)
    
    def to_input_string(self) -> str:
        """
        Converts the CleavageSite object to a string representation for input.

        Returns
        -------
        str
            A string representation of the CleavageSite object in the format:
            C:R1-R2
        """
        return self.Adapter(**(self.model_dump())).to_string()

class CleavageSiteList(BaseObjList):
    
    def describe(self):
        return f"CleavageSiteList with {len(self)} cleavage sites"
    
    def _validate_item(self, item: CleavageSite) -> None:
        """
        Validate that the item is an instance of CleavageSite.
        
        Parameters
        ----------
        item : CleavageSite
            The item to validate.
        
        Raises
        ------
        TypeError
            If the item is not an instance of CleavageSite.
        """
        if not isinstance(item, CleavageSite):
            raise TypeError(f"Item must be an instance of CleavageSite, got {type(item)}")