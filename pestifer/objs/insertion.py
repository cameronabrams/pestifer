# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Insertions are user-specified modifications in which a sequence of one or more amino acids is added to a protein chain.
"""
import logging
logger=logging.getLogger(__name__)

from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.stringthings import split_ri

class Insertion(AncestorAwareObj):
    """
    A class for handling insertions of amino acid residues within an otherwise
    fully resolved chain.
    """

    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum','insertion','sequence']
    """
    Required attributes for an Insertion object.
    These attributes must be provided when creating an Insertion object.
    
    - ``chainID``: The chain ID of the segment where the insertion occurs.
    - ``resseqnum``: The residue number where the insertion is made.
    - ``insertion``: The insertion code for the residue.
    - ``sequence``: The one-letter amino acid sequence to be inserted.
    """
    
    yaml_header='insertions'
    """
    YAML header for Insertion objects.
    This header is used to identify Insertion objects in YAML files.
    """
    
    objcat='seq'
    """
    Category of the Insertion object.
    This categorization is used to group Insertion objects in the object manager.
    """
    
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
