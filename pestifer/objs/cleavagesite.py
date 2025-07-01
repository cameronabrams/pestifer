# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
CleavageSite is a class for handling chain cleavage in molecular structures.
It represents a user-specified modification that cleaves a segment of a chain
between two specified residues, creating a new segment."""
import logging
logger=logging.getLogger(__name__)

from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.stringthings import split_ri

class CleavageSite(AncestorAwareObj):
    """
    A class for handling chain cleavage.  Note that this mod is not expected to be part of an 
    ObjManager so the yaml_header and objcat attributes are irrelevant.  It is instead handled
    as a run task.
    """

    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2']
    """
    Required attributes for a CleavageSite object.
    These attributes must be provided when creating a CleavageSite object.

    - ``chainID``: The chain ID of the segment to be cleaved.
    - ``resseqnum1``: The N-terminal residue number of the cleavage site.
    - ``insertion1``: The insertion code of the N-terminal residue.
    - ``resseqnum2``: The C-terminal residue number of the cleavage site.
    - ``insertion2``: The insertion code of the C-terminal residue.
    """
    
    yaml_header='cleavages'
    """
    YAML header for CleavageSite objects.
    This header is used to identify CleavageSite objects in YAML files.
    """

    objcat='seq'
    """
    Category of the CleavageSite object.
    This categorization is used to group CleavageSite objects in the object manager.
    """

    @singledispatchmethod
    def __init__(self,input_obj):
        super().__init__(input_obj)
    @__init__.register(str)
    def _from_shortcode(self,shortcode):
        # C:R1-R2
        # C: chain ID
        # R1: resid+insertioncode of N-terminal partner in peptide bond
        # R2: resid+insertioncode of C-terminal partner in peptide bond
        items=shortcode.split(':')
        assert len(items)==2,f'Bad cleavage site shortcode: {shortcode}'
        chainID=items[0]
        resrange=items[1]
        ri1,ri2=resrange.split('-')
        r1,i1=split_ri(ri1)
        r2,i2=split_ri(ri2)
        input_dict={
            'chainID':chainID,
            'resseqnum1':r1,
            'insertion1':i1,
            'resseqnum2':r2,
            'insertion2':i2
        }
        super().__init__(input_dict)

class CleavageSiteList(AncestorAwareObjList):
    pass