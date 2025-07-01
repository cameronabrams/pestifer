# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A Substitution is a user-specified modification in which a sequence of one or
more contiguous residues are replaced by a new sequence of one or more residues.
"""
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.stringthings import split_ri

class Substitution(AncestorAwareObj):
    """
    A class for handling substitutions 
    """

    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2','subseq']
    """
    Required attributes for a Substitution object.
    These attributes must be provided when creating a Substitution object.
    
    - ``chainID``: The chain ID of the segment where the substitution occurs.
    - ``resseqnum1``: The N-terminal residue number of the substitution.
    - ``insertion1``: The insertion code of the N-terminal residue.
    - ``resseqnum2``: The C-terminal residue number of the substitution.
    - ``insertion2``: The insertion code of the C-terminal residue.
    - ``subseq``: The one-letter amino acid sequence to be substituted.
    """
    
    yaml_header='substitutions'
    """
    YAML header for Substitution objects.
    This header is used to identify Substitution objects in YAML files.
    """
    
    objcat='seq'
    """
    Category of the Substitution object.
    This categorization is used to group Substitution objects in the object manager.
    """
    
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