# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A Substitution is a user-specified modification in which a sequence of one or
more contiguous residues are replaced by a new sequence of one or more residues.
"""
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from pydantic import Field
from typing import ClassVar, Any

from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.stringthings import split_ri, join_ri

class Substitution(BaseObj):
    """
    A class for handling substitutions 
    """

    _required_fields = {'chainID','resseqnum1','insertion1','resseqnum2','insertion2','subseq'}
    """
    Required attributes for a Substitution object.
    These attributes must be provided when creating a Substitution object.
    
    - ``chainID``: The chain ID of the segment where the substitution occurs.
    - ``resseqnum1``: The N-terminal residue number of the substitution.
    - ``insertion1``: The insertion code of the N-terminal residue.
    - ``resseqnum2``: The C-terminal residue number of the substitution.
    - ``insertion2``: The insertion code of the C-terminal residue.
    - ``subseq``: The one-letter amino acid sequence to be substituted.
    """
    
    chainID: str = Field(..., description="Chain ID of the segment where the substitution occurs")
    resseqnum1: int = Field(..., description="N-terminal residue number of the substitution")
    insertion1: str = Field(..., description="Insertion code of the N-terminal residue")
    resseqnum2: int = Field(..., description="C-terminal residue number of the substitution")
    insertion2: str = Field(..., description="Insertion code of the C-terminal residue")
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
    
    def describe(self):
        return f"Substitution(chainID={self.chainID}, resseqnum1={self.resseqnum1}, insertion1={self.insertion1}, resseqnum2={self.resseqnum2}, insertion2={self.insertion2}, subseq={self.subseq})"

    class Adapter:
        def __init__(self, chainID: str, resseqnum1: int, insertion1: str, resseqnum2: int, insertion2: str, subseq: str):
            self.chainID = chainID
            self.resseqnum1 = resseqnum1
            self.insertion1 = insertion1
            self.resseqnum2 = resseqnum2
            self.insertion2 = insertion2
            self.subseq = subseq

        @classmethod
        def from_string(cls, raw: str):
            # C:nnn-ccc,S
            # where:
            # C is the chain ID
            # nnn is the N-terminal residue/insertion of the sequence to be substituted
            # ccc is the C-terminal residue/insertion of the sequence to be substituted
            # S is the one-letter amino acid sequence to be substituted
            p1=raw.split(':')
            chainID=p1[0]
            p2=p1[1].split(',')
            seqrange=p2[0]
            subseq=p2[1]
            seq=seqrange.split('-')
            r1,i1=split_ri(seq[0])
            r2,i2=split_ri(seq[1])
            input_dict={
                'chainID':chainID,
                'resseqnum1':int(r1),
                'insertion1':i1,
                'resseqnum2':int(r2),
                'insertion2':i2,
                'subseq':subseq
            }
            return cls(**input_dict)
        
        def to_string(self) -> str:
            return f"{self.chainID}:{join_ri(self.resseqnum1, self.insertion1)}-{join_ri(self.resseqnum2, self.insertion2)},{self.subseq}"

        def to_dict(self) -> dict:
            return {k: v for k, v in self.__dict__.items() if not v is None}
    
    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        input_dict = adapter.to_dict()
        return cls(**input_dict)
    
    @singledispatchmethod
    @classmethod
    def new(cls, raw: Any) -> "Substitution":
        raise TypeError(f"Cannot create Substitution from {type(raw)}: {raw}")

    @new.register(str)
    @classmethod
    def _from_string(cls, raw: str) -> "Substitution":
        """
        Create a new Substitution instance from a shortcode string.
        
        Parameters
        ----------
        raw : str
            The shortcode string in the format C_RRR-D_SSS.
        
        Returns
        -------
        Substitution
            A new instance of Substitution.
        """
        adapter = cls.Adapter.from_string(raw)
        instance = cls._from_adapter(adapter)
        return instance
    
    def to_input_string(self) -> str:
        """
        Converts the Substitution object to a string representation for input.
        
        Returns
        -------
        str
            A string representation of the Substitution object in the format:
            C:nnn-ccc,S
        """
        return self.Adapter(**(self.model_dump())).to_string()

class SubstitutionList(BaseObjList[Substitution]):
    """
    A list of Substitution objects.
    This class is used to manage a collection of Substitution objects.
    """
    
    def describe(self):
        return f'<SubstitutionList: {len(self)} items>'
