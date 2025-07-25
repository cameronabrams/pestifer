# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A patch is a modification to a residue or residues defined in the CHARMM force field by ``PRES`` records.
It is used to apply specific modifications to residues in a molecular structure, such as adding or removing
functional groups or modifying the residue's properties.
"""
from copy import deepcopy
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from typing import ClassVar, Optional, Any, List
from pydantic import Field
from ..core.stringthings import split_ri, join_ri
from ..core.baseobj_new import BaseObj, BaseObjList

class Patch(BaseObj):
    """
    A class for handing patch residues
    """
    
    _required_fields = ['patchname', 'chainID', 'resseqnum', 'insertion']
    _optional_fields = ['residue', 'use_in_segment', 'use_after_regenerate']
    """
    Required attributes for a Patch object.
    These attributes must be provided when creating a Patch object.
    
    - ``patchname``: The name of the patch as defined in the CHARMM force field.
    - ``chainID``: The chain identifier of the residue to be patched.
    - ``resseqnum``: The residue sequence number of the residue to be patched.
    - ``insertion``: The insertion code of the residue to be patched.

    Optional attributes for a Patch object.
    - ``residue``: The Residue object to be patched.  Assigned by the PatchList.assign_residues method.
    - ``use_in_segment``: Specifies whether the patch is applied to the first or last residue in a segment, and therefore should be declared inside a psfgen ``segment``.
    - ``use_after_regenerate``: A boolean indicating whether the patch should be applied after the ``regenerate angles dihedrals`` psfgen command.
    """
    patchname: str = Field(..., description="Name of the patch as defined in the CHARMM force field")
    chainID: str = Field(..., description="Chain identifier of the residue to be patched")
    resseqnum: int = Field(..., description="Residue sequence number of the residue to be patched")
    insertion: str = Field(..., description="Insertion code of the residue to be patched")
    residue: Optional[Any] = Field(None, description="Residue object to be patched")
    use_in_segment: Optional[str] = Field(None, description="Specifies whether the patch is applied to the first or last residue in a segment, and therefore should be declared inside a psfgen `segment`")
    use_after_regenerate: Optional[bool] = Field(False, description="If True, the patch is applied after the `regenerate angles dihedrals` psfgen command")

    _yaml_header: ClassVar[str] = 'patches'
    """
    YAML header for Patch objects.
    This header is used to identify Patch objects in YAML files.
    """

    _objcat: ClassVar[str] = 'seq'
    """
    Category of the Patch object.
    This categorization is used to group Patch objects in the object manager.
    """


    def describe(self):
        """
        Describe the Patch object.
        
        Returns
        -------
        str
            A string description of the Patch object, including patch name, chain ID, residue sequence number, insertion code, and whether it is used in a segment or after regeneration.
        """
        return f"Patch(patchname={self.patchname}, chainID={self.chainID}, resseqnum={self.resseqnum}, insertion={self.insertion}, residue={self.residue}, use_in_segment={self.use_in_segment}, use_after_regenerate={self.use_after_regenerate})"

    class Adapter:
        _in_segment_Npatches: ClassVar[List[str]] = ['NTER','GLYP','PROP','ACE','ACED','ACP','ACPD','NNEU','NGNE']
        """
        List of patch names that are applied to the first residue in a segment.
        These patches are typically used for N-terminal modifications.
        """

        _in_segment_Cpatches: ClassVar[List[str]] = ['CTER','CNEU','PCTE','CT1','CT2','CT3']
        """
        List of patch names that are applied to the last residue in a segment.
        These patches are typically used for C-terminal modifications.
        """

        _after_regenerate_patches: ClassVar[List[str]] = ['PHEM','FHEM']
        """
        List of patch names that should be applied after the ``regenerate angles dihedrals`` psfgen command.
        """

        def __init__(self, patchname: str, chainID: str, resseqnum: int, insertion: str, residue: Optional[Any] = None, use_in_segment: Optional[str] = None, use_after_regenerate: bool = False):
            self.patchname = patchname
            self.chainID = chainID
            self.resseqnum = resseqnum
            self.insertion = insertion
            self.residue = residue
            self.use_in_segment = use_in_segment
            self.use_after_regenerate = use_after_regenerate

        @classmethod
        def from_string(cls, raw: str):
            ### shortcode format: patchname:chainID:rrr
            ### patchname: name of the patch in CHARMMFF
            ### chainID: chain identifier of residue to be patched
            ### rrr: residue sequence number of residue to be patched
            parts=raw.split(':')
            if len(parts)<3:
                raise ValueError(f'Invalid patch shortcode: {raw}')
            ri=parts[2]
            r,i=split_ri(ri)
            if parts[0] in cls._in_segment_Npatches:
                use_in_segment='first'
            elif parts[0] in cls._in_segment_Cpatches:
                use_in_segment='last'
            else:
                use_in_segment=''
            if parts[0] in cls._after_regenerate_patches:
                use_after_regenerate=True
            else:
                use_after_regenerate=False
            input_dict={
                'patchname':parts[0],
                'chainID':parts[1],
                'resseqnum':r,
                'insertion':i,
                'residue':None,
                'use_in_segment':use_in_segment,
                'use_after_regenerate':use_after_regenerate
            }
            return cls(**input_dict)
        
        def to_string(self) -> str:
            """
            Convert the Adapter instance to a shortcode string.
            
            Returns
            -------
            str
                The shortcode string in the format patchname:chainID:resseqnum,insertion.
            """
            rescode = join_ri(self.resseqnum, self.insertion)
            return f"{self.patchname}:{self.chainID}:{rescode}"

    @BaseObj.from_input.register(Adapter)
    @classmethod
    def _from_adapter(cls, adapter: Adapter):
        """
        Create a Patch object from an Adapter instance, registered by BaseObj.from_input.
        
        Parameters
        ----------
        adapter : Adapter
            The Adapter instance containing the patch attributes.

        Returns
        -------
        Patch
            A new Patch instance created from the Adapter.
        """
        return cls(
            patchname=adapter.patchname,
            chainID=adapter.chainID,
            resseqnum=adapter.resseqnum,
            insertion=adapter.insertion,
            residue=adapter.residue,
            use_in_segment=adapter.use_in_segment,
            use_after_regenerate=adapter.use_after_regenerate
        )

    @singledispatchmethod
    @classmethod
    def new(cls, raw: Any) -> "Patch":
        pass

    @new.register(str)
    @classmethod
    def _from_str(cls, raw: str) -> "Patch":
        """
        Create a new Patch instance from a shortcode string.

        Parameters
        ----------
        raw : str
            The shortcode string in the format patchname:chainID:resseqnum,insertion.

        Returns
        -------
        Patch
            A new instance of Patch.
        """
        adapter = cls.Adapter.from_string(raw)
        return cls._from_adapter(adapter)

    def to_input_string(self) -> str:
        """
        Convert the Patch object to a string representation for input.
        
        Returns
        -------
        str
            A string representation of the Patch object in the format:
            patchname:chainID:resseqnum,insertion
        """
        return self.Adapter(**(self.model_dump())).to_string()

    def return_TcL(self):
        """
        Writes the patch command.
        """
        return f'patch {self.patchname} {self.chainID}:{self.resseqnum}{self.insertion}'
        
class PatchList(BaseObjList):
    """
    A class for handling lists of Patch objects
    
    This class inherits from BaseObjList and provides methods to manage
    a list of Patch objects.
    """

    def describe(self):
        """
        Returns a string description of the PatchList.
        
        Returns
        -------
        str
            A description of the PatchList, including the number of patches it contains.
        """
        return f'PatchList with {len(self)} patches'
    
    def _validate_item(self, item: Patch) -> None:
        """
        Validate that the item is an instance of Patch.
        
        Parameters
        ----------
        item : Patch
            The item to validate.
        
        Raises
        ------
        TypeError
            If the item is not an instance of Patch.
        """
        if not isinstance(item, Patch):
            raise TypeError(f'Item must be an instance of Patch, got {type(item)} instead.')

    def assign_residues(self, Residues):
        """
        Assigns a list of Residue objects to the patch residues.
        
        Parameters
        ----------
        Residues : list of Residue
            A list of Residue objects to be assigned to the patch residues.
        
        Returns
        -------
        PatchList
            A new PatchList object containing the residues that were not assigned because no applicable residues were found
        """
        logger.debug(f'Patches: Assigning residues from list of {len(Residues)} residues')
        ignored=self.assign_objs_to_attr('residue',Residues,resseqnum='resseqnum',chainID='chainID',insertion='insertion')
        return self.__class__(ignored)
