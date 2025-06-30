# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for processing biological assemblies
"""

import numpy as np
from mmcif.api.PdbxContainers import DataContainer
from pidibble.pdbparse import get_symm_ops
import logging
logger=logging.getLogger(__name__)

from .asymmetricunit import AsymmetricUnit
from ..core.baseobj import AncestorAwareObj, AncestorAwareObjList
from .chainidmanager import ChainIDManager

def build_tmat(RotMat,TransVec):
    """
    Builds a 4 x 4 homogeneous transformation matrix 
    
    Parameters
    ----------
    RotMat: numpy.ndarray
        3 x 3 rotation matrix
    TransVec: numpy.ndarray
        translation vector
    """
    tmat=np.identity(4,dtype=float)
    for i in range(3):
        for j in range(3):
            tmat[i][j]=RotMat[i][j]
        tmat[i][3]=TransVec[i]
    return tmat

class Transform(AncestorAwareObj):
    """
    A class for handling transformations of segments in a molecular structure.
    This class represents a transformation that can be applied to segments
    in an asymmetric unit, including rotation and translation.
    """

    req_attr=AncestorAwareObj.req_attr+['index','tmat','applies_chainIDs','chainIDmap','segname_by_type_map']
    """
    Required attributes for the Transform class.
    
    Attributes
    ----------
    index : int
        Index of the transformation.
    tmat : numpy.ndarray
        4 x 4 transformation matrix.
    applies_chainIDs : list
        List of chain IDs to which this transformation applies.
    chainIDmap : dict
        Mapping of chain IDs for the transformation.
    segname_by_type_map : dict
        Mapping of segment names by type for the transformation.
    """
    def __init__(self,*input_objs):
        input_dict={
            'tmat':None,
            'applies_chainIDs':[],
            'index':'UNSET'
        }
        if len(input_objs)==2: # assume we have a barec from a pdbfile
            barec=input_objs[0]
            input_dict['index']=input_objs[1]
            RotMat,TransVec=get_symm_ops(barec)
            input_dict['applies_chainIDs']=barec.header
            input_dict['tmat']=build_tmat(RotMat,TransVec)
        elif len(input_objs)==4: # M, V, chains, idx
            RotMat=input_objs[0]
            TransVec=input_objs[1]
            input_dict['applies_chainIDs']=input_objs[2]
            input_dict['index']=input_objs[3]
            input_dict['tmat']=build_tmat(RotMat,TransVec)
        elif len(input_objs)>0 and input_objs[0]=='identity':
            RotMat=np.identity(3)
            TransVec=np.zeros(3)
            input_dict['tmat']=build_tmat(RotMat,TransVec)
            input_dict['index']=0
            input_dict['applies_chainIDs']=[]
        input_dict['chainIDmap']={}
        input_dict['segname_by_type_map']={}
        super().__init__(input_dict)

    def is_identity(self):
        """
        Checks if the transformation is an identity transformation.
        An identity transformation is one where the transformation matrix is
        an identity matrix and the translation vector is zero.
        
        Returns
        -------
        bool
            True if the transformation is an identity transformation, False otherwise.
        """
        return np.array_equal(np.identity(4,dtype=float),self.tmat)

    def register_mapping(self,segtype,chainID,seglabel):
        """
        Registers a mapping of segment type to chain ID and segment label.
        This method updates the `segname_by_type_map` attribute with the provided
        segment type, chain ID, and segment label.

        Parameters
        ----------
        segtype : str
            The type of the segment (e.g., 'protein', 'nucleic').
        chainID : str
            The chain ID associated with the segment.
        seglabel : str
            The label for the segment.
        """
        if not segtype in self.segname_by_type_map:
            self.segname_by_type_map[segtype]={}
        self.segname_by_type_map[segtype][chainID]=seglabel

    def write_TcL(self):
        """
        Generates a Tcl command string that represents the transformation.
        This method constructs a string that can be used in a Tcl script to apply
        the transformation to a segment in VMD.

        Returns
        -------
        str
            A string containing the Tcl command to apply the transformation.
        """
        retstr=r'{ '
        for i in range(4):
            retstr+=r'{ '
            for j in range(4):
               retstr+='{} '.format(self.tmat[i][j])
            retstr+=r' } '
        retstr+=r' }'
        return retstr
    
    def __eq__(self,other):
        """
        Checks if this Transform instance is equal to another Transform instance.
        Two Transform instances are considered equal if their transformation matrices are the same.
        
        Parameters
        ----------
        other : Transform
            The other Transform instance to compare with.

        Returns
        -------
        bool
            True if the transformation matrices are equal, False otherwise.
        """
        return np.array_equal(self.tmat,other.tmat)
    
    def generate_chainIDmap(self,auChainIDs,daughters,CM):
        """
        Generates a mapping of chain IDs for the transformation.
        This method creates a mapping of chain IDs that this transformation applies to,
        based on the asymmetric unit's chain IDs and the daughters of segments.
        If the transformation is an identity transformation, it applies a "thru map" to the chain IDs.
        Otherwise, it generates a new mapping for the chain IDs.
        
        Parameters
        ----------
        auChainIDs : list
            The list of chain IDs in the asymmetric unit.
        daughters : dict
            A dictionary mapping parent segment IDs to their daughter segment IDs.
        CM : ChainIDManager
            The ChainIDManager instance used to manage chain ID mappings.
        """
        applies_to=self.applies_chainIDs[:]
        for d,v in daughters.items():
            if d in self.applies_chainIDs:
                applies_to.extend(v)
        logger.debug(f'Transform applies to {applies_to}')
        if self.is_identity():
            logger.debug(f'Identity transform gets a thru map applied to {applies_to}')
            self.chainIDmap=CM.thru_map(auChainIDs,applies_to)
        else:
            logger.debug(f'Transform gets a new map applied to {applies_to}')
            self.chainIDmap=CM.generate_next_map(auChainIDs,applies_to)

class TransformList(AncestorAwareObjList):
    """
    A class for handling lists of Transform objects.
    This class inherits from AncestorAwareObjList and provides methods to manage
    collections of Transform instances.
    """
    def __init__(self,*args):
        L=[]
        if len(args)==1:
            if type(args[0])==AncestorAwareObjList:
                for idx,item in enumerate(args[0]):
                    L.append(Transform(item,idx))
            elif args[0]=='identity':
                L=[Transform('identity')]
        super().__init__(L)

class BioAssemb(AncestorAwareObj):
    """
    A class for handling biological assemblies in molecular structures.
    This class represents a biological assembly, which can consist of multiple
    transformations applied to segments in an asymmetric unit.
    It is initialized with a dictionary or an AsymmetricUnit, TransformList, or AncestorAwareObjList instance.
    If initialized with an AsymmetricUnit, it creates a default assembly with
    the name 'A.U.' and an identity transformation.
    If initialized with a TransformList or AncestorAwareObjList, it uses the transforms from that list.
    If initialized with a dictionary, it expects the dictionary to contain the keys 'name', 'transforms', and 'index'.
    The class also maintains a static index to ensure unique assembly names.
    """
    _index=1 # start at 1
    req_attr=AncestorAwareObj.req_attr+['name','transforms','index']
    """
    Required attributes for the BioAssemb class.
    
    Attributes
    ----------
    name : str
        The name of the biological assembly.
    transforms : TransformList
        A list of Transform objects representing the transformations applied to the assembly.
    index : int
        An index for the biological assembly, used to ensure unique names.
    """

    def __init__(self,input_obj):
        if type(input_obj)==dict:
            input_dict=input_obj
        elif type(input_obj)==AsymmetricUnit:
            input_dict={
                'name':'A.U.',
                'transforms':TransformList('identity')
            }
        elif type(input_obj)==AncestorAwareObjList:
            input_dict={'transforms':TransformList(input_obj)}
        elif type(input_obj)==TransformList:
            input_dict={'transforms':input_obj}
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
        """
        Activate the biological assembly by generating chain ID maps for its transformations.
        
        Parameters
        ----------
        AU : AsymmetricUnit
            The asymmetric unit to which this biological assembly applies.
        CM : ChainIDManager
            The ChainIDManager instance used to manage chain ID mappings.
        """
        for T in self.transforms:
            T.generate_chainIDmap(AU.segments.segnames,AU.segments.daughters,CM)

class BioAssembList(AncestorAwareObjList):
    """
    A class for handling lists of BioAssemb objects.
    This class inherits from AncestorAwareObjList and provides methods to manage
    collections of biological assemblies.
    """

    def __init__(self,*obj):
        BioAssemb.reset_index()
        B=[]
        if len(obj)==1:
            p_struct=obj[0] # whole darn thing
            if type(p_struct)==dict:
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
                    reclist=AncestorAwareObjList([])
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
                    B.append(BioAssemb(reclist))
                logger.debug(f'There are {len(B)} biological assemblies')
            elif type(p_struct)==DataContainer:
                Assemblies=p_struct.getObj('pdbx_struct_assembly')
                gen=p_struct.getObj('pdbx_struct_assembly_gen')
                oper=p_struct.getObj('pdbx_struct_oper_list')
                for ba_idx in range(len(Assemblies)):
                    logger.debug(f'CIF: Establishing BA {ba_idx}')
                    assemb_id=Assemblies.getValue('id',ba_idx)
                    this_gen_idx_list=gen.selectIndices(assemb_id,'assembly_id')
                    logger.debug(f'BA {ba_idx} points to {len(this_gen_idx_list)} gen indexes')
                    transforms=TransformList()
                    for this_gen_idx in this_gen_idx_list:
                        this_oper_list=gen.getValue('oper_expression',this_gen_idx).split(',')
                        logger.debug(f'BA {ba_idx} gen {this_gen_idx} opers {this_oper_list}')
                        this_asyms=gen.getValue('asym_id_list',this_gen_idx).split(',')
                        logger.debug(f'asym ids: {this_asyms}')
                        idx=0
                        # logger.debug(f'Expecting {len(this_opers)} transforms')
                        for k,opere in enumerate(this_oper_list):
                            oper_idx=oper.selectIndices(opere,'id')[0]
                            logger.debug(f'making transform from oper {oper_idx}')
                            m=np.identity(3)
                            v=np.zeros(3)
                            for i in range(3):
                                I=i+1
                                vlabel=f'vector[{I}]'
                                v[i]=float(oper.getValue(vlabel,oper_idx))
                                for j in range(3):
                                    J=j+1
                                    mlabel=f'matrix[{I}][{J}]'
                                    m[i][j]=float(oper.getValue(mlabel,oper_idx))
                            T=Transform(m,v,this_asyms,idx)
                            transforms.append(T)
                            idx+=1
                    logger.debug(f'parsed {len(transforms)} transforms for ba {ba_idx}')
                    BA=BioAssemb(transforms)
                    B.append(BA)
                logger.debug(f'There '+'is' if len(B)==1 else 'are'+f' {len(B)} biological assembl'+'y' if len(B)==1 else 'ies')
        super().__init__(B)


