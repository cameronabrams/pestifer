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
    _required_fields = {'chainID', 'resseqnum1', 'insertion1', 'resseqnum2', 'insertion2'}

    chainID:    str = Field(..., description="Chain ID of the segment to be cleaved")
    resseqnum1: int = Field(..., description="N-terminal residue number of the cleavage site")
    insertion1: str = Field(..., description="Insertion code of the N-terminal residue")
    resseqnum2: int = Field(..., description="C-terminal residue number of the cleavage site")
    insertion2: str = Field(..., description="Insertion code of the C-terminal residue")

    _yaml_header: ClassVar[str] = 'cleavages'
    _objcat: ClassVar[str] = 'seq'

    @staticmethod
    def _adapt(*args) -> dict:
        if isinstance(args[0], str):
            return CleavageSite._from_shortcode(args[0])
        raise TypeError(f"Cannot convert {type(args[0])} to CleavageSite")

    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        chainID, res_range = raw.split(":")
        resseqnum1, insertion1 = split_ri(res_range.split("-")[0])
        resseqnum2, insertion2 = split_ri(res_range.split("-")[1])
        return dict(
            chainID=chainID,
            resseqnum1=resseqnum1,
            insertion1=insertion1,
            resseqnum2=resseqnum2,
            insertion2=insertion2
        )

    def shortcode(self) -> str:
        return f"{self.chainID}:{join_ri(self.resseqnum1, self.insertion1)}-{join_ri(self.resseqnum2, self.insertion2)}"

class CleavageSiteList(BaseObjList[CleavageSite]):
    def describe(self):
        return f"CleavageSiteList with {len(self)} cleavage sites"
