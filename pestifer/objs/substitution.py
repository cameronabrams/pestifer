# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A Substitution is a user-specified modification in which a sequence of one or
more contiguous residues are replaced by a new sequence of one or more residues.
"""
import logging
logger = logging.getLogger(__name__)

from pydantic import Field
from typing import ClassVar

from ..core.baseobj import BaseObj, BaseObjList
from .resid import ResID, ResIDList

class Substitution(BaseObj):
    """
    A class for handling substitutions 
    """

    _required_fields = {'chainID','resid1','resid2','subseq'}
    """
    Required attributes for a Substitution object.
    These attributes must be provided when creating a Substitution object.
    
    - ``chainID``: The chain ID of the segment where the substitution occurs.
    - ``resid1``: The N-terminal residue ID of the substitution.
    - ``resid2``: The C-terminal residue ID of the substitution.
    - ``subseq``: The one-letter amino acid sequence to be substituted.
    """
    
    chainID: str = Field(..., description="Chain ID of the segment where the substitution occurs")
    resid1: ResID = Field(..., description="N-terminal residue ID of the substitution")
    resid2: ResID = Field(..., description="C-terminal residue ID of the substitution")
    subseq: str = Field(..., description="One-letter amino acid sequence to be substituted")

    _yaml_header: ClassVar[str] ='substitutions'
    """
    YAML header for Substitution objects.
    This header is used to identify Substitution objects in YAML files.
    """

    _objcat: ClassVar[str] = 'seq'
    """
    Category of the Substitution object.
    This categorization is used to group Substitution objects in the object manager.
    """
    
    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for Substitution instantiation.
        This method is used to convert various input formats into a dictionary that can be used to create a Substitution object.
        """
        if args and len(args) == 1 and isinstance(args[0], str):
            # C:nnn-ccc,S
            # where:
            # C is the chain ID
            # nnn is the N-terminal residue/insertion of the sequence to be substituted
            # ccc is the C-terminal residue/insertion of the sequence to be substituted
            # S is the one-letter amino acid sequence to be substituted
            raw = args[0]
            if ':' not in raw or '-' not in raw or ',' not in raw:
                raise ValueError(f'Invalid substitution shortcode: {raw}')
            p1 = raw.split(':')
            chainID = p1[0]
            p2 = p1[1].split(',')
            seqrange = p2[0]
            subseq = p2[1]
            seq = ResIDList(seqrange)
            resid1, resid2 = seq
            input_dict = {
                'chainID': chainID,
                'resid1': resid1,
                'resid2': resid2,
                'subseq': subseq
            }
            return input_dict
        return super()._adapt(*args, **kwargs)

    def shortcode(self) -> str:
        return f"{self.chainID}:{self.resid1.resid}-{self.resid2.resid},{self.subseq}"

class SubstitutionList(BaseObjList[Substitution]):
    """
    A list of Substitution objects.
    This class is used to manage a collection of Substitution objects.
    """
    
    def describe(self):
        return f'<SubstitutionList: {len(self)} items>'
