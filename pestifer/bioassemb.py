

import numpy as np
from pidibble.pdbparse import get_symm_ops
import logging
logger=logging.getLogger(__name__)
from .basemod import AncestorAwareMod, AncestorAwareModList

class BiomT(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['index','tmat']
    @classmethod
    def from_rot_trans(cls,RotMat:np.ndarray,TransVec:np.ndarray,index):
        tmat=np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0]],dtype=float)
        for i in range(3):
            for j in range(3):
                tmat[i][j]=RotMat[i][j]
            tmat[i][3]=TransVec[i]
        input_dict={
            'tmat':tmat,
            'index':index
        }
        return cls(input_dict)
    def write_TcL(self):
        retstr=r'{ '
        for i in range(3):
            retstr+=r'{ '
            for j in range(4):
               retstr+='{} '.format(self.tmat[i][j])
            retstr+=r' } '
        retstr+='{ 0 0 0 1 } }'
        return retstr
    def __eq__(self,other):
        return np.array_equal(self.tmat,other.tmat)

class BiomTList(AncestorAwareModList):
    @classmethod
    def from_pdbrecord(cls,pdbrecord):
        L=[]
        M,T=get_symm_ops(pdbrecord)
        for i,(m,t) in enumerate(zip(M,T)):
            L.append(BiomT.from_rot_trans(m,t,i))
        return cls(L)

class BioAssemb(AncestorAwareMod):
    _index=1 # start at 1
    req_attr=AncestorAwareMod.req_attr+['name','chainIDs','biomt','index']
    ''' Container for handling info for "REMARK 350 BIOMOLECULE: #" stanzas in RCSB PDB files
        or _pdbx_struct blocks in mmCIF files '''
    def __init__(self,input_dict):
        input_dict['index']=BioAssemb._index
        BioAssemb._index+=1
        super().__init__(input_dict)

    @classmethod
    def from_pdbrecord(cls,pdbrecord):
        rs=pdbrecord.key.split('.')
        input_dict={
            'name':'.'.join(rs[1:3]),
            'chainIDs':pdbrecord.header,
            'biomt':BiomTList.from_pdbrecord(pdbrecord) 
        }
        inst=cls(input_dict)
        return inst

class BioAssembList(AncestorAwareModList):
    @classmethod
    def from_pdb(cls,p_struct):
        B=[]
        barecs=[p_struct[x] for x in p_struct if ('REMARK.350.BIOMOLECULE' in x and 'TRANSFORM' in x)]
        for rec in barecs:
            B.append(BioAssemb.from_pdbrecord(rec))
        return cls(B)
