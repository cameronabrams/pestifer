# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod
from ..stringthings import split_ri

from ..baseobj import AncestorAwareObj, AncestorAwareObjList

class Patch(AncestorAwareObj):
    """A class for handing patch residues
    
    Attributes
    ----------
    req_attr: list
        * patchname: (str) name of the patch in CHARMMFF
        * chainID: (str) chain identifier of residue to be patched
        * resseqnum: (int) residue number of residue to be patched
        * insertion: (chr) residue insertion code of residue to be patched
    
    yaml_header: (str)
        label used for yaml format input/output
    
    objcat: (str)
        type of mod from among objcats

    """
    req_attr=['patchname','chainID','resseqnum','insertion','residue']
    yaml_header='patches'
    objcat='seq'

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
        input_dict={
            'patchname':parts[0],
            'chainID':parts[1],
            'resseqnum':r,
            'insertion':i,
            'residue':None
        }
        super().__init__(input_dict)
        
class PatchList(AncestorAwareObjList):
 def assign_residues(self,Residues):
        logger.debug(f'Patches: Assigning residues from list of {len(Residues)} residues')
        ignored=self.assign_objs_to_attr('residue',Residues,resseqnum='resseqnum',chainID='chainID',insertion='insertion')
        return self.__class__(ignored)
