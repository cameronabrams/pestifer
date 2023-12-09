#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Bonds
"""
from .mods import *
from pidibble.baserecord import BaseRecord
from functools import singledispatchmethod

class Bond(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['atom_idx1','atom_idx2']
    opt_attr=AncestorAwareMod.opt_attr+['atom1','atom2']
    yaml_header='atoms'
    PDB_keyword='ATOM'
    mmCIF_name='atom_site'

class BondList(AncestorAwareModList):
    @singledispatchmethod
    def __init__(self,inputobj):
        super().__init__(inputobj)
    
    @__init__.register(str)
    def _from_psffile(self,filename):
        pass

    @__init__.register(dict)
    def _from_bonddict(self,bonddict):
        pass