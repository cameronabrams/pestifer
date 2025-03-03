# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)

from functools import singledispatchmethod

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..stringthings import split_ri

class Insertion(AncestorAwareObj):
    """A class for handling insertions of amino acid residues within an otherwise
    fully resolved chain
    
    Attributes
    ----------
    req_attr: list
        * chainID: chain id where insertion is to be made
        * resseqnum, insertion: resid+ins.code marking position AFTER which inserted residues are to be inserted
        * sequence: sequence of one-letter amino acids defining the insertion
    """
    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum','insertion','sequence']
    yaml_header='insertions'
    objcat='seq'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # C,R,SSS
        # C: chain ID
        # R: resid+insertioncode
        # SSS: one-letter amino acid sequence to insert
        items=shortcode.split(',')
        assert len(items)==3,f'Bad insertion shortcode: {shortcode}'
        rescode=items[1]
        self.integer_increment=False
        if rescode[-1]=='+':
            self.integer_increment=True
            rescode=rescode[:-1]
        r,i=split_ri(rescode)
        input_dict={
            'chainID':items[0],
            'resseqnum':r,
            'insertion':i,
            'sequence':items[2]
        }
        super().__init__(input_dict)

class InsertionList(AncestorAwareObjList):
    pass
