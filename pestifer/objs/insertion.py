# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Insertions are user-specified modifications in which a sequence of one or more amino acids is added to a protein chain.
"""
import logging
logger=logging.getLogger(__name__)

from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.stringthings import split_ri, join_ri
from typing import Optional, ClassVar
from pydantic import Field

class Insertion(BaseObj):
    """
    A class for handling insertions of amino acid residues within an otherwise
    fully resolved chain.
    """
    _required_fields = ['chainID', 'resseqnum', 'insertion', 'sequence']
    _optional_fields = ['integer_increment']

    chainID: str = Field(..., description="Chain ID of the segment where the insertion occurs")
    resseqnum: int = Field(..., description="Residue sequence number where the insertion is made")
    insertion: str = Field(..., description="Insertion code for the residue")
    sequence: str = Field(..., description="One-letter amino acid sequence to be inserted")
    integer_increment: Optional[bool] = Field(False, description="If True, the residue sequence number is incremented by 1 for each insertion; indicated by a '+' at the end of the residue code in the shortcode format.")
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

    def describe(self):
        """
        Describe the Insertion object.
        
        Returns
        -------
        str
            A string description of the Insertion object, including chain ID, residue number, insertion code, and sequence.
        """
        return f"Insertion(chainID={self.chainID}, resseqnum={self.resseqnum}, insertion={self.insertion}, sequence={self.sequence}, integer_increment={self.integer_increment})"
    
    class Adapter:
        """
        A class to represent the shortcode format for Insertion, so that we can register to BaseObj.from_input rather than defining a local from_input.
        
        The shortcode format is C,R,SSS
        where:
        - C is the chain ID
        - R is the residue number and insertion code
        - SSS is the one-letter amino acid sequence to insert after R
        """
        def __init__(self, chainID: str, resseqnum: int, insertion: str, sequence: str, integer_increment: bool = False):
            self.chainID = chainID
            self.resseqnum = resseqnum
            self.insertion = insertion
            self.sequence = sequence
            self.integer_increment = integer_increment

        @classmethod
        def from_string(cls, raw: str):
            items = raw.split(',')
            if len(items) != 3:
                raise ValueError(f"Bad insertion shortcode: {raw}")
            rescode = items[1]
            integer_increment = False
            if rescode[-1] == '+':
                integer_increment = True
                rescode = rescode[:-1]
            r, i = split_ri(rescode)
            return cls(items[0], r, i, items[2], integer_increment)

        def to_string(self) -> str:
            rescode = join_ri(self.resseqnum, self.insertion)
            if self.integer_increment:
                rescode += '+'
            return f"{self.chainID},{rescode},{self.sequence}"
        
    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        """
        Create an Insertion object from an adapter.
        
        Parameters
        ----------
        adapter : Adapter
            An instance of Insertion.Adapter containing the necessary attributes to create an Insertion object.
        
        Returns
        -------
        Insertion
            An instance of Insertion created from the adapter.
        """
        return cls(
            chainID=adapter.chainID,
            resseqnum=adapter.resseqnum,
            insertion=adapter.insertion,
            sequence=adapter.sequence,
            integer_increment=adapter.integer_increment
        )
    
    @classmethod
    def new(cls, raw: str) -> "Insertion":
        """
        Create a new Insertion instance from a shortcode string.
        
        Parameters
        ----------
        raw : str
            The shortcode string in the format C,R,SSS.
        
        Returns
        -------
        Insertion
            A new Insertion instance.
        """
        adapter = cls.Adapter.from_string(raw)
        return cls._from_adapter(adapter)

class InsertionList(BaseObjList[Insertion]):
    def describe(self) -> str:
        """
        Describe the InsertionList.
        
        Returns
        -------
        str
            A string description of the InsertionList, including the number of insertions.
        """
        return f"InsertionList with {len(self)} insertions"
    
    def _validate_item(self, item: Insertion) -> None:
        """
        Validate that the item is an instance of Insertion.
        
        Parameters
        ----------
        item : Insertion
            The item to validate.
        
        Raises
        ------
        TypeError
            If the item is not an instance of Insertion.
        """
        if not isinstance(item, Insertion):
            raise TypeError(f"Item must be an instance of Insertion, got {type(item)}")
