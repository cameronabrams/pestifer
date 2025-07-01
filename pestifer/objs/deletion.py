# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A deletion is a user-specified modification in which a sequence of one or
more contiguous residues are deleted.
"""
import logging
logger=logging.getLogger(__name__)
from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.stringthings import split_ri

class Deletion(AncestorAwareObj):
    """
    A class for handling deletions in a molecular sequence.
    """

    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2']
    """
    Required attributes for a Deletion object.
    These attributes must be provided when creating a Deletion object.

    - ``chainID``: The chain ID of the segment from which residues are deleted.
    - ``resseqnum1``: The N-terminal residue number of the deletion.
    - ``insertion1``: The insertion code of the N-terminal residue.
    - ``resseqnum2``: The C-terminal residue number of the deletion.
    - ``insertion2``: The insertion code of the C-terminal residue.
    """

    yaml_header='deletions'
    """
    YAML header for Deletion objects.
    This header is used to identify Deletion objects in YAML files.
    """
    
    objcat='seq'
    """
    Category of the Deletion object.
    This categorization is used to group Deletion objects in the object manager.
    """
    
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
