#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Some classes/methods for interfacing with mmCIF files.
"""
import logging

from pidibble.pdbparse import PDBParser
from mmcif.io.IoAdapterCore import IoAdapterCore
from mmcif.api.PdbxContainers import DataContainer
from collections import UserDict
from pathlib import Path

logger=logging.getLogger(__name__)

class CIFdict(UserDict):
    """
    A class for generating a custom-format dictionary from an mmcif input object
    """
    def __init__(self,Obj,idx,lowercase=True,blankers=[' ','?']):
        if lowercase:
            data={c.lower():Obj.getValue(c,idx) for c in Obj.getAttributeList()}
        else:
            data={c:Obj.getValue(c,idx) for c in Obj.getAttributeList()}
        for k,v in data.items():
            if v in blankers:
                data[k]=''
        super().__init__(data)

def CIFload(input_obj) -> DataContainer:
    """
    Downloads (if necessary) and reads in a mmCIF file into an mmcif DataContainer object
    """
    filepath=input_obj
    if type(input_obj) is str and '.pdb' not in input_obj:
        filepath=f'{input_obj}.cif'
        # fetch the cif file if not already present
        PDBParser(PDBcode=input_obj,input_format='mmCIF').fetch()
    elif type(input_obj) is Path:
        PDBParser(filepath=input_obj, input_format='mmCIF').fetch()
        filepath=input_obj
    # strip out the pdbx_audit_revision_item that confuses
    # vmd
    logger.debug(f'Stripping the offensive \'pdbx_audit_revision_item\' block from {filepath} and resaving')
    dataList=IoAdapterCore().readFile(filepath)
    dataList[0].remove('pdbx_audit_revision_item')
    IoAdapterCore().writeFile(filepath, dataList)
    # return the first DataContainer in the list that .readFile() returns
    return IoAdapterCore().readFile(filepath)[0]

