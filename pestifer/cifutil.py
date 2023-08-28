
class CIFobject:
    pass

class CIFdict(dict):
    def __init__(self,Obj,idx):
        self.data={c:Obj.getValue(c,idx) for c in Obj.getAttributeList()}

    def __getitem__(self,key):
        return self.data[key]
        # this will make bool give false if data is empty...

    def __bool__(self):
        return True
    
    def get(self,key,default):
        if key in self.data:
            return self.data[key]
        else:
            return default


def CIFMakeStructs(db):
    structs={}
    for k in db.keys():
        kk=k.split('.')
        s=kk[0]
        if s in structs:
            i+=1
            structs[s][kk[1]]=i
        else:
            i=0
            structs[s]={}
            structs[s][kk[1]]=i
    return structs

def GetCIFStructDictList(db,structs,struk):
    keys=[_ for _ in structs[struk].keys()]
    dlist=[]
    chek=struk+'.'+keys[0]
    isloop=True if db.FindLoop(chek)!=-1 else False
    if isloop:
        x=db.GetLoop(chek)
        for y in x:
            d={}
            for v in keys:
                d[v]=y[structs[struk][v]]
            dlist.append(d)
    else:
        d={}
        for v in keys:
            d[v]=db[struk+'.'+v]
        dlist.append(d)
    return dlist