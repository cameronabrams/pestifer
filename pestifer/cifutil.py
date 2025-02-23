#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Some classes/methods for interfacing with mmcif
"""
from pidibble.pdbparse import PDBParser
from mmcif.io.IoAdapterCore import IoAdapterCore
from collections import UserDict

class CIFdict(UserDict):
    """A class for generating a custom-format dictionary from an mmcif input object"""
    def __init__(self,Obj,idx,lowercase=True,linkers={'.':{'label_seq_id':'auth_seq_id'}},blankers=[' ','?']):
        if lowercase:
            data={c.lower():Obj.getValue(c,idx) for c in Obj.getAttributeList()}
        else:
            data={c:Obj.getValue(c,idx) for c in Obj.getAttributeList()}
        for k,v in data.items():
            klkey=linkers.get(v,{}).get(k,k)
            valk=data[klkey]
            if valk in blankers:
                valk=''
            data[klkey]=valk
        super().__init__(data)

def CIFload(pdb_id):
    """Downloads (if necessary) and reads in a mmCIF file into an mmcif DataContainer object"""
    # fetch the cif file if not already present
    PDBParser(PDBcode=pdb_id,input_format='mmCIF').fetch()
    # return the first DataContainer in the list that .readFile() returns
    return IoAdapterCore().readFile(f'{pdb_id}.cif')[0]
