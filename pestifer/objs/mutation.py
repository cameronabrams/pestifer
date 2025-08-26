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

from pydantic import Field
from typing import ClassVar

from .resid import ResID
from .seqadv import Seqadv
from ..core.baseobj import BaseObj, BaseObjList
from ..core.labels import Labels

logger = logging.getLogger(__name__)

class Mutation(BaseObj):
    """
    A class for handling single-residue mutations
    """

    _required_fields = {'chainID', 'origresname', 'resid', 'newresname', 'typekey'}
    """
    Required attributes for a Mutation object.
    These attributes must be provided when creating a Mutation object.

    - ``chainID``: The chain ID of the segment where the mutation occurs.
    - ``origresname``: The original residue name at the mutation site.
    - ``resid``: The residue ID where the mutation is made, represented by a `ResID` object.
    - ``newresname``: The new residue name to be mutated to.
    - ``typekey``: A key indicating the type of mutation (e.g., ``user``, ``author``, ``mmCIF``).
    """

    _optional_fields = {'pdbx_auth_seq_num'}
    """
    Optional attributes for a Mutation object.
    These attributes may be present but are not required.

    - ``pdbx_auth_seq_num``: The PDBx author sequence number, which is used in mmCIF files to indicate the residue number.
    """

    chainID: str = Field(..., description="Chain ID of the segment where the mutation occurs")
    origresname: str = Field(..., description="Original residue name at the mutation site")
    resid: ResID = Field(..., description="Residue ID where the mutation is made")
    newresname: str = Field(..., description="New residue name to be mutated to")
    typekey: str = Field(..., description="Key indicating the type of mutation (e.g., 'user', 'author', 'mmCIF')")
    pdbx_auth_seq_num: int | None = Field(None, description="PDBx author sequence number, used in mmCIF files")

    _yaml_header: ClassVar[str] = 'mutations'
    """
    YAML header for Mutation objects.
    This header is used to identify Mutation objects in YAML files.
    """

    _objcat: ClassVar[str] = 'seq'
    """
    Category of the Mutation object.
    This categorization is used to group Mutation objects in the object manager.
    """

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for Mutation instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if args and isinstance(args[0], str):
            input_dict = Mutation._from_shortcode(args[0])
            return input_dict
        elif args and isinstance(args[0], Seqadv):
            input_dict = Mutation._from_seqadv(args[0])
            return input_dict
        return super()._adapt(*args, **kwargs)
    
    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        """
        Converts a shortcode string to a dictionary of parameters for Mutation.
        The shortcode format is: C:nnn,rrr,mmm
        where:
        - C is the chain ID
        - nnn is the original residue name
        - rrr is the residue sequence number plus optional insertion code
        - mmm is the new residue name
        """
        s1 = raw.split(':')
        chainID = s1[0]
        mut_spec = s1[1]
        # possible formats for mut_spec
        #
        # 1. nnn,rrr,mmm
        # 2. nnnrrrmmm
        # 3. n,rrr,m
        # 4. nrrrm
        #
        # nnn, mmm are three-letter residue codes
        # rrr is a resid with possible insertion code
        # n, m are one-letter residue codes
        if ',' not in mut_spec:
            # we have a comma-less string
            # if the string is using one-letter resid codes, then the second
            # character will be a digit
            if mut_spec[1].isdigit():
                idx = 1
            else:
                idx = 3
                # assume the first letter is a residue code
            origresname = Labels.res_123.get(mut_spec[:idx],mut_spec[:idx])
            # assume the last letter is a residue code
            newresname = Labels.res_123.get(mut_spec[-idx:],mut_spec[-idx:])
            # assume letters/numbers in the middle are a resid in string form
            resid = ResID(mut_spec[idx:-idx])
        else:
            s2 = s1[1].split(',')
            origresname = Labels.res_123.get(s2[0], s2[0])
            resid = ResID(s2[1])
            newresname = Labels.res_123.get(s2[2], s2[2])
        typekey = 'user'  # Default type key
        return dict(
            chainID=chainID,
            origresname=origresname,
            resid=resid,
            newresname=newresname,
            typekey=typekey
        )
    
    @staticmethod
    def _from_seqadv(seqadv: Seqadv) -> dict:
        """
        Converts a Seqadv object to a dictionary of parameters for Mutation.
        
        Parameters
        ----------
        seqadv : Seqadv
            The Seqadv object containing mutation information.
        
        Returns
        -------
        dict
            A dictionary with the attributes needed to create a Mutation object.
        """
        return dict(
            chainID=seqadv.chainID,
            origresname=seqadv.resname,
            resid=ResID(seqadv.resid),
            newresname=seqadv.dbRes,
            typekey=seqadv.typekey,
            pdbx_auth_seq_num=seqadv.pdbx_auth_seq_num
        )
    
    def shortcode(self) -> str:
        """
        Returns a shortcode representation of the Mutation object.

        Returns
        -------
        str
            A string in the format C:nnn,rrr,mmm where C is the chain ID, nnn is the original residue name,
            rrr is the residue sequence number plus optional insertion code, and mmm is the new residue name.
        """
        return f"{self.chainID}:{self.origresname},{self.resid.resid},{self.newresname}"

    def write_TcL(self):
        """
        Writes the Tcl command to perform the mutation in a Psfgen script.
        
        Returns
        -------
        str
            The Tcl command to mutate the residue at the specified position to the new residue name.
        """
        # if hasattr(self,'pdbx_auth_seq_num'): # mmCIF, no insertion code
        #     return f'    mutate {self.resid.resseqnum} {self.newresname}'
        # else:
        return f'    mutate {self.resid.resid} {self.newresname}'

class MutationList(BaseObjList[Mutation]):
    def describe(self):
        return f'<MutationList with {len(self)} mutations>'
