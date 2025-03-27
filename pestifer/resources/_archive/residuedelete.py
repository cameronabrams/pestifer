# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ...baseobj import AncestorAwareObj, AncestorAwareObjList
from ...stringthings import split_ri

class ResidueDelete(AncestorAwareObj):
    """A class for handling residue deletions
    
    Attributes
    ----------
    req_attr: list
        * chainID: chain identifier of residue to be deleted
        * resseqnum: resid of residue to be deleted
        * insertion: insertion code of residue to be deleted

    """
    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum','insertion']
    yaml_header='residuedeletions'
    objcat='seq'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: C_nnn
        # C -- chainID
        # nnn -- resid of residue to be deleted
        p1=shortcode.split('_')
        chainID=p1[0]
        ri=p1[1]
        r,i=split_ri(ri)
        input_dict={
            'chainID':chainID,
            'resseqnum':int(r),
            'insertion':i
        }
        super().__init__(input_dict)

class ResidueDeleteList(AncestorAwareObjList):
    pass