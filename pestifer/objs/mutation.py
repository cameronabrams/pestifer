# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Single-residue mutations are handled in psfgen by including a "mutate" directive within
a protein segment.  The mutate directive needs only the resid and 3-letter desired
residue name at that resid position.  Mutate directives come after all pdb and residue
directives in segment, typically.

Mutations can be inferred directly from coordinate-file metadata or supplied by the
user.
"""
import logging
import pdb
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from typing import ClassVar, Any
from pydantic import Field
from .seqadv import Seqadv

from ..core.baseobj_new import BaseObj, BaseObjList
from ..core.labels import Labels
from ..core.stringthings import split_ri, join_ri

class Mutation(BaseObj):
    """
    A class for handling single-residue mutations
    """

    _required_fields = ['chainID', 'origresname', 'resseqnum', 'insertion', 'newresname', 'typekey']
    """
    Required attributes for a Mutation object.
    These attributes must be provided when creating a Mutation object.

    - ``chainID``: The chain ID of the segment where the mutation occurs.
    - ``origresname``: The original residue name at the mutation site.
    - ``resseqnum``: The residue number where the mutation is made.
    - ``insertion``: The insertion code for the original residue.
    - ``newresname``: The new residue name to be mutated to.
    - ``typekey``: A key indicating the type of mutation (e.g., ``user``, ``author``, ``mmCIF``).
    """

    _optional_fields = ['pdbx_auth_seq_num']
    """
    Optional attributes for a Mutation object.
    These attributes may be present but are not required.

    - ``pdbx_auth_seq_num``: The PDBx author sequence number, which is used in mmCIF files to indicate the residue number.
    """

    chainID: str = Field(..., description="Chain ID of the segment where the mutation occurs")
    origresname: str = Field(..., description="Original residue name at the mutation site")
    resseqnum: int = Field(..., description="Residue sequence number where the mutation is made")
    insertion: str = Field(..., description="Insertion code for the original residue")
    newresname: str = Field(..., description="New residue name to be mutated to")
    typekey: str = Field(..., description="Key indicating the type of mutation (e.g., 'user', 'author', 'mmCIF')")
    pdbx_auth_seq_num: int = Field(None, description="PDBx author sequence number, used in mmCIF files")

    _yaml_header: ClassVar[str] =   'mutations'
    """
    YAML header for Mutation objects.
    This header is used to identify Mutation objects in YAML files.
    """

    _objcat: ClassVar[str] = 'seq'
    """
    Category of the Mutation object.
    This categorization is used to group Mutation objects in the object manager.
    """

    def describe(self):
        """
        Describe the Mutation object.
        
        Returns
        -------
        str
            A string description of the Mutation object, including chain ID, original residue name, residue sequence number, insertion code, new residue name, and type key.
        """
        return f"Mutation(chainID={self.chainID}, origresname={self.origresname}, resseqnum={self.resseqnum}, insertion={self.insertion}, newresname={self.newresname}, typekey={self.typekey}, pdbx_auth_seq_num={self.pdbx_auth_seq_num})"
    
    class Adapter:
        def __init__(self, chainID: str, origresname: str, resseqnum: int, insertion: str, newresname: str, typekey: str):
            self.chainID = chainID
            self.origresname = origresname
            self.resseqnum = resseqnum
            self.insertion = insertion
            self.newresname = newresname
            self.typekey = typekey
        
        @classmethod
        def from_string(cls, raw: str):
            """
            Create an Adapter instance from a shortcode string.
            
            Parameters
            ----------
            raw : str
                The shortcode string in the format C:nnn,rrr,mmm.
                where C is the chain ID, nnn is the original residue name, rrr is the residue sequence number plus optional insertion code, mmm is the new residue name.

            Returns
            -------
            Adapter
                An Adapter instance initialized with the parsed values from the shortcode string.
            """
            s1 = raw.split(':')
            chainID = s1[0]
            s2 = s1[1].split(',')
            origresname = s2[0]
            resseqnum,insertion = split_ri(s2[1])
            newresname = s2[2]
            typekey = 'user'  # Default type key

            return cls(chainID=chainID, origresname=origresname, resseqnum=resseqnum, insertion=insertion, newresname=newresname, typekey=typekey)

        @classmethod
        def from_seqadv(cls, seqadv: Seqadv):
            """
            Create an Adapter instance from a Seqadv object.

            Parameters
            ----------
            seqadv : Seqadv
                The Seqadv object containing mutation information.

            Returns
            -------
            Adapter
                An Adapter instance initialized with the values from the Seqadv object.
            """
            return cls(chainID=seqadv.chainID, origresname=seqadv.resname, resseqnum=seqadv.resseqnum, insertion=seqadv.insertion, newresname=seqadv.dbRes, typekey=seqadv.typekey)


    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        """
        Create a Mutation object from an Adapter instance, registered by BaseObj.from_input.
        
        Parameters
        ----------
        adapter : Adapter
            An instance of Mutation.Adapter containing the necessary attributes to create a Mutation object.
        
        Returns
        -------
        Mutation
            A new Mutation instance initialized with the attributes from the adapter.
        """
        input_dict = adapter.__dict__
        return cls(**input_dict)

    @singledispatchmethod
    @classmethod
    def new(cls, raw: Any) -> "Mutation":
        pass

    @new.register(str)
    @classmethod
    def _from_shortcode(cls, raw: str) -> "Mutation":
        return cls._from_adapter(cls.Adapter.from_string(raw))
    
    @new.register(Seqadv)
    @classmethod
    def _from_seqadv(cls, seqadv: Seqadv) -> "Mutation":
        """
        Create a new Mutation instance from a Seqadv object.
        
        Parameters
        ----------
        seqadv : Seqadv
            The Seqadv object containing mutation information.
        
        Returns
        -------
        Mutation
            A new Mutation instance initialized with the values from the Seqadv object.
        """
        adapter = cls.Adapter.from_seqadv(seqadv)
        return cls._from_adapter(adapter)

    def __str__(self):
        return f'{self.chainID}:{self.origresname}{self.resseqnum}{self.insertion}{self.newresname}'

    def write_TcL(self):
        """
        Writes the Tcl command to perform the mutation in a Psfgen script.
        
        Returns
        -------
        str
            The Tcl command to mutate the residue at the specified position to the new residue name.
        """
        if hasattr(self,'pdbx_auth_seq_num'): # mmCIF!
            return f'    mutate {self.resseqnum} {self.newresname}'
        else:
            resseqnumi=f'{self.resseqnum}{self.insertion}'
            return f'    mutate {resseqnumi} {self.newresname}'

class MutationList(BaseObjList[Mutation]):
    def describe(self):
        return f'MutationList with {len(self)} mutations'
    
    def _validate_item(self, item:Mutation):
        if not isinstance(item, Mutation):
            raise TypeError(f"Expected Mutation, got {type(item)}")
        
