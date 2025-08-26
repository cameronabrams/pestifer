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

from ..core.baseobj import BaseObj, BaseObjList
from .resid import ResID

class CleavageSite(BaseObj):
    """
    A class for handling chain cleavage.  Note that this mod is not expected to be part of an 
    ObjManager so the yaml_header and objcat attributes are irrelevant.  It is instead handled
    as a run task.
    """
    _required_fields = {'chainID', 'resid1', 'resid2'}

    chainID: str = Field(..., description="Chain ID of the segment to be cleaved")
    resid1: ResID = Field(..., description="N-terminal residue information of the cleavage site")
    resid2: ResID = Field(..., description="C-terminal residue information of the cleavage site")

    _yaml_header: ClassVar[str] = 'cleavages'
    _objcat: ClassVar[str] = 'seq'

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        if args and isinstance(args[0], str):
            return CleavageSite._from_shortcode(args[0])
        return super()._adapt(*args, **kwargs)

    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        chainID, res_range = raw.split(":")
        resid1 = ResID(res_range.split("-")[0])
        resid2 = ResID(res_range.split("-")[1])
        return dict(
            chainID=chainID,
            resid1=resid1,
            resid2=resid2
        )

    def shortcode(self) -> str:
        return f"{self.chainID}:{self.resid1.resid}-{self.resid2.resid}"

class CleavageSiteList(BaseObjList[CleavageSite]):
    def describe(self):
        return f"CleavageSiteList with {len(self)} cleavage sites"
