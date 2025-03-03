# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..stringthings import split_ri

class Deletion(AncestorAwareObj):
    """A class for handling deletions 
    
    A deletion is a user-specified modification in which a sequence of one or
    more contiguous residues are deleted.

    Attributes
    ----------
    req_attr: list
        * chainID: chain identifier
        * resseqnum1: N-terminal resid of sequence to be deleted
        * insertion1: insertion code of N-terminal residue
        * resseqnum2: C-terminal resid of sequence to be deleted
        * insertion2: insertion code of the C-terminal residue

    """
    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2']
    # opt_attr=AncestorAwareObj.opt_attr+['model']
    yaml_header='deletions'
    objcat='seq'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: C:nnn-ccc
        # C -- chainID
        # nnn -- N-terminal resid of sequence to be deleted
        # ccc -- C-terminal resid of sequence to be deleted
        p1=shortcode.split(':')
        chainID=p1[0]
        seq=p1[1].split('-')
        r1,i1=split_ri(seq[0])
        r2,i2=split_ri(seq[1])          
        input_dict={
            'chainID':chainID,
            'resseqnum1':int(r1),
            'insertion1':i1,
            'resseqnum2':int(r2),
            'insertion2':i2
        }
        super().__init__(input_dict)

class DeletionList(AncestorAwareObjList):
    pass
