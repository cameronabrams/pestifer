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
from .resid import ResID

class Deletion(BaseObj):
    """
    A class for handling deletions in a molecular sequence.
    """
    _required_fields = {'chainID', 'resid1', 'resid2'}

    chainID: str = Field(..., description="Chain ID of the segment from which residues are deleted")
    resid1: ResID = Field(..., description="N-terminal residue of the deletion")
    resid2: ResID = Field(..., description="C-terminal residue of the deletion")

    _yaml_header: ClassVar[str] = 'deletions'
    _objcat: ClassVar[str] = 'seq'
    
    @staticmethod
    def _adapt(*args) -> dict:
        """
        Adapts the input to a dictionary format suitable for Deletion instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if isinstance(args[0], str):
            input_dict = Deletion._from_shortcode(args[0])
            return input_dict
        raise TypeError(f"Cannot convert {type(args[0])} to Deletion")

    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        """
        Converts a shortcode string to a dictionary of parameters for Deletion.
        The shortcode format is: C:nnn-ccc
        where:
        - C is the chain ID
        - nnn is the N-terminal residue number of the deletion
        - ccc is the C-terminal residue number of the deletion
        """
        chainID, res_range = raw.split(":")
        resid1 = ResID(res_range.split("-")[0])
        resid2 = ResID(res_range.split("-")[1])
        return dict(
            chainID = chainID,
            resid1 = resid1,
            resid2 = resid2
        )

    def shortcode(self) -> str:
        return f"{self.chainID}:{self.resid1.resid}-{self.resid2.resid}"

class DeletionList(BaseObjList[Deletion]):
    def describe(self):
        return f'DeletionList with {len(self)} deletions'
