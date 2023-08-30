from pidibble.pdbparse import PDBParser
from mmcif.io.IoAdapterCore import IoAdapterCore
from collections import UserDict

class CIFdict(UserDict):
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

# def CIFMakeStructs(db):
#     structs={}
#     for k in db.keys():
#         kk=k.split('.')
#         s=kk[0]
#         if s in structs:
#             i+=1
#             structs[s][kk[1]]=i
#         else:
#             i=0
#             structs[s]={}
#             structs[s][kk[1]]=i
#     return structs

# def GetCIFStructDictList(db,structs,struk):
#     keys=[_ for _ in structs[struk].keys()]
#     dlist=[]
#     chek=struk+'.'+keys[0]
#     isloop=True if db.FindLoop(chek)!=-1 else False
#     if isloop:
#         x=db.GetLoop(chek)
#         for y in x:
#             d={}
#             for v in keys:
#                 d[v]=y[structs[struk][v]]
#             dlist.append(d)
#     else:
#         d={}
#         for v in keys:
#             d[v]=db[struk+'.'+v]
#         dlist.append(d)
#     return dlist