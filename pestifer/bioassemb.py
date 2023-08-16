

import numpy as np
from pidibble.pdbparse import get_symm_ops, PDBRecord
from pidibble.baserecord import BaseRecord
import logging
logger=logging.getLogger(__name__)
from .basemod import AncestorAwareMod, AncestorAwareModList
from .asymmetricunit import AsymmetricUnit

class BiomT(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['index','tmat','chainIDmap','segname_by_type_map']
    
    def __init__(self,*input_objs):
        if len(input_objs)==3:
            RotMat,TransVec,index=input_objs
        else:
            RotMat=np.identity(3)
            TransVec=np.zeros(3)
            index=0
        tmat=np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0]],dtype=float)
        for i in range(3):
            for j in range(3):
                tmat[i][j]=RotMat[i][j]
            tmat[i][3]=TransVec[i]
        input_dict={
            'tmat':tmat,
            'index':index
        }
        input_dict['chainIDmap']={}
        input_dict['segname_by_type_map']={}
        super().__init__(input_dict)

    def register_mapping(self,segtype,chainID,seglabel):
        if not segtype in self.segname_by_type_map:
            self.segname_by_type_map[segtype]={}
        self.segname_by_type_map[segtype][chainID]=seglabel

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
    def __init__(self,*args):
        L=[]
        if len(args)==1:
            if type(args[0])==PDBRecord or type(args[0])==BaseRecord:
                pdbrecord=args[0]
                M,T=get_symm_ops(pdbrecord)
                L=[BiomT(m,t,i) for i,(m,t) in enumerate(zip(M,T))]
        else:
            L=[BiomT()]
        super().__init__(L)

class BioAssemb(AncestorAwareMod):
    _index=1 # start at 1
    req_attr=AncestorAwareMod.req_attr+['name','chainIDs','biomt','index']
    ''' Container for handling info for "REMARK 350 BIOMOLECULE: #" stanzas in RCSB PDB files '''
    def __init__(self,input_obj):
        if type(input_obj)==dict:
            input_dict=input_obj
        elif type(input_obj) in [PDBRecord,BaseRecord]:
            pdbrecord=input_obj
            rs=pdbrecord.key.split('.')
            input_dict={
                'name':'.'.join(rs[1:3]),
                'chainIDs':pdbrecord.header,
                'biomt':BiomTList(pdbrecord) 
            }
        elif type(input_obj)==AsymmetricUnit:
            au=input_obj
            input_dict={
                'name':'A.U.',
                'chainIDs':au.chainIDs,
                'biomt':BiomTList()
            }
        else:
            logger.warning(f'Cannot initialize {type(self)} from object of type {type(input_obj)}')
        input_dict['index']=BioAssemb._index
        BioAssemb._index+=1
        super().__init__(input_dict)
    @classmethod
    def reset_index(cls):
        cls._index=1

class BioAssembList(AncestorAwareModList):
    def __init__(self,*obj):
        BioAssemb.reset_index()
        B=[]
        if len(obj)==1:
            p_struct=obj[0]
            barecs=[p_struct[x] for x in p_struct if ('REMARK.350.BIOMOLECULE' in x and 'TRANSFORM' in x)]
            B=[BioAssemb(rec) for rec in barecs]
        super().__init__(B)
