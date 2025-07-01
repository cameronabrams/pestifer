# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A patch is a modification to a residue or residues defined in the CHARMM force field by ``PRES`` records.
It is used to apply specific modifications to residues in a molecular structure, such as adding or removing
functional groups or modifying the residue's properties.
"""
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from ..core.stringthings import split_ri
from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList

class Patch(AncestorAwareObj):
    """
    A class for handing patch residues
    """
    
    req_attr=['patchname','chainID','resseqnum','insertion','residue','use_in_segment']
    """
    Required attributes for a Patch object.
    These attributes must be provided when creating a Patch object.
    
    - ``patchname``: The name of the patch as defined in the CHARMM force field.
    - ``chainID``: The chain identifier of the residue to be patched.
    - ``resseqnum``: The residue sequence number of the residue to be patched.
    - ``insertion``: The insertion code of the residue to be patched.
    - ``residue``: The Residue object to be patched.
    - ``use_in_segment``: Specifies whether the patch is applied to the first or last residue in a segment, and therefore should be declared inside a psfgen ``segment``.
    """

    yaml_header='patches'
    """
    YAML header for Patch objects.
    This header is used to identify Patch objects in YAML files.
    """
    
    objcat='seq'
    """
    Category of the Patch object.
    This categorization is used to group Patch objects in the object manager.
    """
    
    in_segment_Npatches=['NTER','GLYP','PROP','ACE','ACED','ACP','ACPD','NNEU','NGNE']
    """
    List of patch names that are applied to the first residue in a segment.
    These patches are typically used for N-terminal modifications.
    """
    
    in_segment_Cpatches=['CTER','CNEU','PCTE','CT1','CT2','CT3']
    """
    List of patch names that are applied to the last residue in a segment.
    These patches are typically used for C-terminal modifications.
    """

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
        if not hasattr(self,'insertion'):
            self.insertion=''

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        ### shortcode format: patchname:chainID:rrr
        ### patchname: name of the patch in CHARMMFF
        ### chainID: chain identifier of residue to be patched
        ### rrr: residue sequence number of residue to be patched
        parts=shortcode.split(':')
        if len(parts)<3:
            raise ValueError(f'Invalid patch shortcode: {shortcode}')
        ri=parts[2]
        r,i=split_ri(ri)
        if parts[0] in self.in_segment_Npatches:
            use_in_segment='first'
        elif parts[0] in self.in_segment_Cpatches:
            use_in_segment='last'
        else:
            use_in_segment=''
        input_dict={
            'patchname':parts[0],
            'chainID':parts[1],
            'resseqnum':r,
            'insertion':i,
            'residue':None,
            'use_in_segment':use_in_segment
        }
        super().__init__(input_dict)
        
class PatchList(AncestorAwareObjList):
    """
    A class for handling lists of Patch objects
    
    This class inherits from AncestorAwareObjList and provides methods to manage
    a list of Patch objects.
    """
    def assign_residues(self,Residues):
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
