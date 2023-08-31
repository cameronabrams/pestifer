from pidibble.pdbparse import PDBParser
from mmcif.io.IoAdapterCore import IoAdapterCore
from collections import UserDict

class CIFdict(UserDict):
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
    PDBParser(PDBcode=pdb_id,input_format='mmCIF').fetch()
    return IoAdapterCore().readFile(f'{pdb_id}.cif')[0]
