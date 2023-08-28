

import numpy as np
from pidibble.pdbparse import get_symm_ops, PDBRecord
from pidibble.baserecord import BaseRecord
import logging
logger=logging.getLogger(__name__)
from .basemod import AncestorAwareMod, AncestorAwareModList
from .asymmetricunit import AsymmetricUnit
from .chainids import ChainIDManager

class Transform(AncestorAwareMod):
    req_attr=AncestorAwareMod.req_attr+['index','tmat','applies_chainIDs','chainIDmap','segname_by_type_map']
    def __init__(self,*input_objs):
        if len(input_objs)==2:
            barec=input_objs[0]
            index=input_objs[1]
            RotMat,TransVec=get_symm_ops(barec)
            applies_chainIDs=barec.header
        else:
            RotMat=np.identity(3)
            TransVec=np.zeros(3)
            index=0
            applies_chainIDs=[]
        tmat=np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0]],dtype=float)
        for i in range(3):
            for j in range(3):
                tmat[i][j]=RotMat[i][j]
            tmat[i][3]=TransVec[i]
        input_dict={
            'tmat':tmat,
            'applies_chainIDs':applies_chainIDs,
            'index':index
        }
        input_dict['chainIDmap']={}
        input_dict['segname_by_type_map']={}
        super().__init__(input_dict)

    def is_identity(self):
        tmat=np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0]],dtype=float)
        return np.array_equal(tmat,self.tmat)

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
    
    def generate_chainIDmap(self,auChainIDs,chainIDmap):
        if self.is_identity():
            logger.debug(f'Identity transform gets a thru map')
            self.chainIDmap=chainIDmap.thru_map(auChainIDs,self.applies_chainIDs)
        else:
            self.chainIDmap=chainIDmap.generate_next_map(auChainIDs,self.applies_chainIDs)

class TransformList(AncestorAwareModList):
    def __init__(self,*args):
        L=[]
        if len(args)==1:
            if type(args[0])==AncestorAwareModList:
                for idx,item in enumerate(args[0]):
                    L.append(Transform(item,idx))
        else:
            L=[Transform()]
        super().__init__(L)

class BioAssemb(AncestorAwareMod):
    _index=1 # start at 1
    # req_attr=AncestorAwareMod.req_attr+['name','chainIDs','biomt','index']
    req_attr=AncestorAwareMod.req_attr+['name','transforms','index']
    ''' Container for handling info for "REMARK 350 BIOMOLECULE: #" stanzas in RCSB PDB files '''
    def __init__(self,input_obj,actual_index=-1):
        if type(input_obj)==dict:
            input_dict=input_obj
        elif type(input_obj)==AsymmetricUnit:
            input_dict={
                'name':'A.U.',
                'transforms':TransformList()
            }
        elif type(input_obj)==AncestorAwareModList:
            input_dict={'transforms':TransformList(input_obj)}
        else:
            logger.warning(f'Cannot initialize {type(self)} from object of type {type(input_obj)}')
        if 'name' not in input_dict:
            input_dict['name']=f'Assembly{BioAssemb._index}'
        input_dict['index']=BioAssemb._index
        BioAssemb._index+=1
        super().__init__(input_dict)

    @classmethod
    def reset_index(cls):
        cls._index=1

    def activate(self,AU:AsymmetricUnit,CM:ChainIDManager):
        for T in self.transforms:
            T.generate_chainIDmap(AU.chainIDs,CM)

class BioAssembList(AncestorAwareModList):
    def __init__(self,*obj):
        BioAssemb.reset_index()
        B=[]
        if len(obj)==1:
            p_struct=obj[0] # whole darn thing
            # extract any ba records
            bareclabels=[x for x in p_struct if ('REMARK.350.BIOMOLECULE' in x and 'TRANSFORM' in x)]
            tr={}
            for lab in bareclabels:
                fs=lab.split('.')
                assert fs[0]=='REMARK'
                assert fs[1]=='350'
                assert 'BIOMOLECULE' in fs[2]
                assert 'TRANSFORM' in fs[3]
                banumber=int(fs[2][11:])
                trnumber=int(fs[3][9:])
                if not banumber in tr:
                    tr[banumber]=[]
                tr[banumber].append(trnumber)
            for ba,trs in tr.items():
                reclist=AncestorAwareModList([])
                savhdr=[]
                for t in trs:
                    if f'REMARK.350.BIOMOLECULE{ba}.TRANSFORM{t}' in p_struct:
                        barec=p_struct[f'REMARK.350.BIOMOLECULE{ba}.TRANSFORM{t}']
                        if hasattr(barec,'header'):
                            savhdr=barec.header
                        else:
                            barec.header=savhdr
                        reclist.append(barec)
                        logger.debug(f'BA {ba} header {barec.header}')
                        logger.debug(barec.pstr())
                B.append(BioAssemb(reclist,actual_index=ba))
            logger.debug(f'There are {len(B)} biological assemblies')
        super().__init__(B)
