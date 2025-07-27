# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A deletion is a user-specified modification in which a sequence of one or
more contiguous residues are deleted.
"""
import logging
logger=logging.getLogger(__name__)

from pydantic import Field
from typing import ClassVar
from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.stringthings import split_ri, join_ri

class Deletion(BaseObj):
    """
    A class for handling deletions in a molecular sequence.
    """
    _required_fields = {'chainID', 'resseqnum1', 'insertion1', 'resseqnum2', 'insertion2'}

    chainID: str = Field(..., description="Chain ID of the segment from which residues are deleted")
    resseqnum1: int = Field(..., description="N-terminal residue number of the deletion")
    insertion1: str = Field(..., description="Insertion code of the N-terminal residue")
    resseqnum2: int = Field(..., description="C-terminal residue number of the deletion")
    insertion2: str = Field(..., description="Insertion code of the C-terminal residue")

    _yaml_header: ClassVar[str] = 'deletions'
    _objcat: ClassVar[str] = 'seq'
    
    class Adapter:
        """
        A class to represent the shortcode format for Deletion, so that we can register to BaseObj.from_input rather than defining a local from_input.

        The shortcode format is C:nnn-ccc
        where:
        - C is the chain ID
        - nnn is the N-terminal residue of the sequence to be deleted
        - ccc is the C-terminal residue of the sequence to be deleted
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
        
        def to_dict(self) -> dict:
            return {k: v for k, v in self.__dict__.items() if not v is None}

    def describe(self):
        return f"Deleteion(chainID={self.chainID}, resseqnum1={self.resseqnum1}, insertion1={self.insertion1}, resseqnum2={self.resseqnum2}, insertion2={self.insertion2})"

    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_shortcode(cls, shortcode: Adapter):
        return cls(**(shortcode.to_dict()))

    @classmethod
    def new(cls, raw: str) -> "Deletion":
        """
        Create a new Deletion instance from a shortcode string.
        
        Parameters
        ----------
        raw : str
            The shortcode string in the format C:nnn-ccc.
        
        Returns
        -------
        Deletion
            A new Deletion instance.
        """
        adapter = cls.Adapter.from_string(raw)
        return cls._from_shortcode(adapter)

    def to_input_string(self) -> str:
        """
        Convert the Deletion object to a string representation for input.
        
        Returns
        -------
        str
            A string representation of the Deletion object in the format:
            C:nnn-ccc
        """
        return self.Adapter(**(self.model_dump())).to_string()

class DeletionList(BaseObjList[Deletion]):
    def describe(self):
        return f'DeletionList with {len(self)} deletions'
