# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
logger=logging.getLogger(__name__)

from functools import singledispatchmethod

from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from ..core.stringthings import split_ri

class CleavageSite(AncestorAwareObj):
    """A class for handling chain cleavage.  Note that this mod is not expected to be part of an 
    ObjManager so the yaml_header and objcat attributes are irrelevant.  It is instead handled
    as a run task.
    
    Attributes
    ----------
    req_attr : list
        * chainID : str
            chain ID of the segment to which cleavage is made
        * resseqnum1 : int
            N-terminal resid of cleavage site
        * insertion1 : str
            insertion code of N-terminal residue
        * resseqnum2 : int
            C-terminal resid of cleavage site
        * insertion2 : str
            insertion code of the C-terminal residue"""
    req_attr=AncestorAwareObj.req_attr+['chainID','resseqnum1','insertion1','resseqnum2','insertion2']
    yaml_header='cleavages'
    objcat='seq'

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