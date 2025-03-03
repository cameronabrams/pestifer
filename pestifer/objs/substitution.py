# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ..baseobj import AncestorAwareObj, AncestorAwareObjList
from ..stringthings import split_ri

class Substitution(AncestorAwareObj):
    """A class for handling substitutions 
    
    A substitution is a user-specified modification in which a sequence of one or
    more contiguous residues are replaced by a new sequence of one or more residues.

    Attributes
    ----------
    req_attr: list
        * chainID: chain identifier
        * resseqnum1: N-terminal resid of sequence to be replaced
        * insertion1: insertion code of N-terminal residue
        * resseqnum2: C-terminal resid of sequence to be replaced
        * insertion2: insertion code of the C-terminal residue
        * subseq: 1-byte rescode sequence to be substituted in

    """
    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2','subseq']
    yaml_header='substitutions'
    objcat='seq'

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)

    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # shortcode format: C:nnn-ccc,abcdef
        # C -- chainID
        # nnn -- N-terminal resid of sequence to be replaced
        # ccc -- C-terminal resid of sequence to be replaced
        # abcdef -- the 1-byte rescode sequence to be substituted
        p1=shortcode.split(':')
        chainID=p1[0]
        p2=p1[1].split(',')
        seqrange=p2[0]
        subseq=p2[1]
        seq=seqrange.split('-')
        r1,i1=split_ri(seq[0])
        r2,i2=split_ri(seq[1])
        input_dict={
            'chainID':chainID,
            'resseqnum1':int(r1),
            'insertion1':i1,
            'resseqnum2':int(r2),
            'insertion2':i2,
            'subseq':subseq
        }
        super().__init__(input_dict)

class SubstitutionList(AncestorAwareObjList):
    pass