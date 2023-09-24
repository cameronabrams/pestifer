#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Some classes/methods for interfacing with mmcif
"""
from pidibble.pdbparse import PDBParser
from mmcif.io.IoAdapterCore import IoAdapterCore
from collections import UserDict

class CIFdict(UserDict):
    """A class for generating a custom-format dictionary from an mmcif input object"""
    # values shown as keys in valuemap are mapped to valuemap values in the resulting dict
    valuemap={'?':'','.':''}
    def __init__(self,Obj,idx,lowercase=True):
        if lowercase:
            data={c.lower():Obj.getValue(c,idx) for c in Obj.getAttributeList()}
        else:
            data={c:Obj.getValue(c,idx) for c in Obj.getAttributeList()}
        for k,v in data.items():
            if v in self.valuemap.keys():
                data[k]=self.valuemap[v]
        super().__init__(data)

def CIFload(pdb_id):
    """Downloads (if necessary) and reads in a mmCIF file into an mmcif object"""
    PDBParser(PDBcode=pdb_id,input_format='mmCIF').fetch()
    return IoAdapterCore().readFile(f'{pdb_id}.cif')[0]
