# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Insertions are user-specified modifications in which a sequence of one or more amino acids is added to a protein chain.
"""
import logging
logger=logging.getLogger(__name__)

from ..core.baseobj_new import BaseObj, BaseObjList
from typing import ClassVar
from pydantic import Field
from .resid import ResID

class Insertion(BaseObj):
    """
    A class for handling insertions of amino acid residues within an otherwise
    fully resolved chain.
    """
    _required_fields = {'chainID', 'resid', 'sequence'}
    _optional_fields = {'integer_increment'}

    chainID: str = Field(..., description="Chain ID of the segment where the insertion occurs")
    resid: ResID = Field(..., description="Residue ID where the insertion is made")
    sequence: str = Field(..., description="One-letter amino acid sequence to be inserted")
    integer_increment: bool | None = Field(None, description="If True, the residue sequence number is incremented by 1 for each insertion; indicated by a '+' at the end of the residue code in the shortcode format.")
    """
    Required attributes for an Insertion object.
    These attributes must be provided when creating an Insertion object.
    
    - ``chainID``: The chain ID of the segment where the insertion occurs.
    - ``resseqnum``: The residue number where the insertion is made.
    - ``insertion``: The insertion code for the residue.
    - ``sequence``: The one-letter amino acid sequence to be inserted.

    Optional attributes for an Insertion object.
    These attributes can be provided to modify the behavior of the insertion.
    - ``integer_increment``: Optional boolean indicating if the residue sequence number is incremented by 1 for each residue substituted.
    """

    _yaml_header: ClassVar[str] = 'insertions'
    """
    YAML header for Insertion objects.
    This header is used to identify Insertion objects in YAML files.
    """

    _objcat: ClassVar[str] = 'seq'
    """
    Category of the Insertion object.
    This categorization is used to group Insertion objects in the object manager.
    """

    @staticmethod
    def _adapt(*args) -> dict:
        """
        Adapts the input to a dictionary format suitable for Insertion instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if isinstance(args[0], str):
            input_dict = Insertion._from_shortcode(args[0])
            return input_dict
        raise TypeError(f"Cannot convert {type(args[0])} to Insertion")

    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        """
        Converts a shortcode string to a dictionary of parameters for Insertion.
        The shortcode format is: C,R,SSS
        where:
        - C is the chain ID
        - R is the residue number and insertion code
        - SSS is the one-letter amino acid sequence to insert after R
        """
        items = raw.split(',')
        if len(items) != 3:
            raise ValueError(f"Bad insertion shortcode: {raw}")
        rescode = items[1]
        integer_increment = False
        if rescode[-1] == '+':
            integer_increment = True
            rescode = rescode[:-1]
        r = ResID(rescode)
        return {
            'chainID': items[0],
            'resid': r,
            'sequence': items[2],
            'integer_increment': integer_increment
        }

    def shortcode(self) -> str:
        """
        Returns the shortcode representation of the Insertion object.
        
        The shortcode format is: C,R,SSS
        where:
        - C is the chain ID
        - R is the residue number and insertion code
        - SSS is the one-letter amino acid sequence to insert after R
        """
        if self.integer_increment:
            rescode += '+'
        return f"{self.chainID},{self.resid.resid},{self.sequence}"

class InsertionList(BaseObjList[Insertion]):
    def describe(self) -> str:
        """
        Describe the InsertionList.
        
        Returns
        -------
        str
            A string description of the InsertionList, including the number of insertions.
        """
        return f"<InsertionList with {len(self)} insertions>"

