# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A patch is a modification to a residue or residues defined in the CHARMM force field by ``PRES`` records.
It is used to apply specific modifications to residues in a molecular structure, such as adding or removing
functional groups or modifying the residue's properties.
"""
from __future__ import annotations

import logging
from typing import ClassVar, TYPE_CHECKING
from pydantic import Field
from ..core.baseobj import BaseObj, BaseObjList
from .resid import ResID

if TYPE_CHECKING:
    from ..molecule.residue import Residue, ResidueList

logger = logging.getLogger(__name__)

class Patch(BaseObj):
    """
    A class for handing patch residues
    """

    _required_fields = {'patchname', 'chainID', 'resid'}
    """
    Required attributes for a Patch object.
    These attributes must be provided when creating a Patch object.
    
    - ``patchname``: The name of the patch as defined in the CHARMM force field.
    - ``chainID``: The chain identifier of the residue to be patched.
    - ``resid``: The residue ID of the residue to be patched.
    """
    _optional_fields = {'residue', 'use_in_segment', 'use_after_regenerate'}
    """
    Optional attributes for a Patch object.
    - ``residue``: The Residue object to be patched.  Assigned by the PatchList.assign_residues method.
    - ``use_in_segment``: Specifies whether the patch is applied to the first or last residue in a segment, and therefore should be declared inside a psfgen ``segment``.
    - ``use_after_regenerate``: A boolean indicating whether the patch should be applied after the ``regenerate angles dihedrals`` psfgen command.
    """
    patchname: str = Field(..., description="Name of the patch as defined in the CHARMM force field")
    chainID: str = Field(..., description="Chain identifier of the residue to be patched")
    resid: ResID = Field(..., description="Residue ID of the residue to be patched")
    residue: "Residue" = Field(None, description="Residue object to be patched")
    use_in_segment: str | None = Field(None, description="Specifies whether the patch is applied to the first or last residue in a segment, and therefore should be declared inside a psfgen `segment`")
    use_after_regenerate: bool | None = Field(False, description="If True, the patch is applied after the `regenerate angles dihedrals` psfgen command")

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

    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        parts = raw.split(':')
        if len(parts) < 3:
            raise ValueError(f'Invalid patch shortcode: {raw}')
        resid = ResID(parts[2])
        if parts[0] in Patch._in_segment_Npatches:
            use_in_segment = 'first'
        elif parts[0] in Patch._in_segment_Cpatches:
            use_in_segment = 'last'
        else:
            use_in_segment = ''
        if parts[0] in Patch._after_regenerate_patches:
            use_after_regenerate = True
        else:
            use_after_regenerate = False
        return {
            'patchname': parts[0],
            'chainID': parts[1],
            'resid': resid,
            'use_in_segment': use_in_segment,
            'use_after_regenerate': use_after_regenerate
        }
    
    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for Patch instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if args and isinstance(args[0], str):
            return Patch._from_shortcode(args[0])
        return super()._adapt(*args, **kwargs)

    _in_segment_Npatches: ClassVar[list[str]] = ['NTER', 'GLYP', 'PROP', 'ACE', 'ACED', 'ACP', 'ACPD', 'NNEU', 'NGNE']
    """
    List of patch names that are applied to the first residue in a segment.
    These patches are typically used for N-terminal modifications.
    """

    _in_segment_Cpatches: ClassVar[list[str]] = ['CTER', 'CNEU', 'PCTE', 'CT1', 'CT2', 'CT3']
    """
    List of patch names that are applied to the last residue in a segment.
    These patches are typically used for C-terminal modifications.
    """

    _after_regenerate_patches: ClassVar[list[str]] = ['PHEM', 'FHEM']
    """
    List of patch names that should be applied after the ``regenerate angles dihedrals`` psfgen command.
    """

    def shortcode(self) -> str:
        """
        Convert the Adapter instance to a shortcode string.
        
        Returns
        -------
        str
            The shortcode string in the format patchname:chainID:resID.
        """
        return f"{self.patchname}:{self.chainID}:{self.resid.resid}"

    def return_TcL(self):
        """
        Writes the patch command.
        """
        return f'patch {self.patchname} {self.chainID}:{self.resid.resid}'

    def assign_residue(self, residue: 'Residue'):
        """
        Assigns a Residue object to the patch.
        """
        self.residue = residue

class PatchList(BaseObjList[Patch]):
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
        return f'<PatchList with {len(self)} patches>'

    def assign_residues(self, Residues: "ResidueList") -> PatchList:
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
        ignored = self.assign_objs_to_attr('residue', Residues, resid='resid', chainID='chainID')
        return self.__class__(ignored)
