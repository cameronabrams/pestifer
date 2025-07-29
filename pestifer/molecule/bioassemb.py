# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for processing biological assemblies
"""

import numpy as np
from functools import singledispatchmethod
from mmcif.api.PdbxContainers import DataContainer
from pidibble.pdbparse import get_symm_ops
from pidibble.pdbrecord import PDBRecord
import logging
logger=logging.getLogger(__name__)
from typing import List, Dict, Optional, Any, ClassVar
from collections import UserList

from .asymmetricunit import AsymmetricUnit
from .chainidmanager import ChainIDManager
from ..util.coord import build_tmat

class Transform:
    """
    A class for handling transformations of segments in a molecular structure.
    This class represents a transformation that can be applied to segments
    in an asymmetric unit, including rotation and translation.

    This method used the :func:`pidibble.pdbparse.get_symm_ops` function to extract the rotation matrix and translation vector from a PDB record.
    It constructs a 4 x 4 transformation matrix from these components.
    
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
    def __init__(self, index: Optional[int] = 0, tmat: Optional[np.ndarray] = None, applies_chainIDs: Optional[List[str]] = None, chainIDmap: Optional[Dict[str, str]] = None, segname_by_type_map: Optional[Dict[str, str]] = None):
        self.index = index
        self.tmat = tmat if tmat is not None else np.identity(4, dtype=float)
        self.applies_chainIDs = applies_chainIDs if applies_chainIDs is not None else []
        self.chainIDmap = chainIDmap if chainIDmap is not None else {}
        self.segname_by_type_map = segname_by_type_map if segname_by_type_map is not None else {}

    @classmethod
    def from_pdb_record(cls, barec: PDBRecord, index: int = 0):
        RotMat,TransVec=get_symm_ops(barec)
        tmat = build_tmat(RotMat, TransVec)
        applies_chainIDs = barec.header if hasattr(barec, 'header') else []
        return cls(index=index, tmat=tmat, applies_chainIDs=applies_chainIDs)
    
    @classmethod
    def from_matvec(cls, RotMat: np.ndarray, TransVec: np.ndarray, applies_chainIDs: List[str] = [], index: int = 0):
        tmat = build_tmat(RotMat, TransVec)
        return cls(index=index, tmat=tmat, applies_chainIDs=applies_chainIDs)

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
        return np.array_equal(np.identity(4,dtype=float), self.tmat)

    def register_mapping(self, segtype, chainID, seglabel):
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
        return np.array_equal(self.tmat, other.tmat)

    def generate_chainIDmap(self, auChainIDs: List[str], daughters: dict, CM: ChainIDManager):
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
            self.chainIDmap=CM.thru_map(auChainIDs, applies_to)
        else:
            logger.debug(f'Transform gets a new map applied to {applies_to}')
            self.chainIDmap=CM.generate_next_map(auChainIDs, applies_to)

class TransformList(UserList):
    """
    A class for handling lists of Transform objects.
    This class inherits from AncestorAwareObjList and provides methods to manage
    collections of Transform instances.
    """
    @classmethod
    def identity(cls):
        """
        Creates an identity TransformList.
        """
        return cls([Transform('identity')])
            
    @classmethod
    def from_pdb_records(cls, pdb_records: List[PDBRecord]):
        """
        Creates a TransformList from a list of PDB records.

        Parameters
        ----------
        pdb_records : list
            A list of PDB records from which to create Transform instances.

        Returns
        -------
        TransformList
            An instance of TransformList containing the created Transform objects.
        """
        transforms = [Transform.from_pdb_record(barec, index=i) for i, barec in enumerate(pdb_records)]
        return cls(transforms)

class BioAssemb:
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

    _index: ClassVar[int] = 1  # start at 1

    @singledispatchmethod
    def __init__(self, input_obj: Any):
        raise TypeError(f'Cannot initialize {type(self)} from object of type {type(input_obj)}')
    
    @__init__.register(dict)
    def _from_dict(self, input_obj: Dict[str, Any]):
        BioAssemb._index+=1
        self.name = input_obj.get('name', f'Assembly{BioAssemb._index}')
        self.transforms = input_obj.get('transforms', TransformList.identity())

    @__init__.register(AsymmetricUnit)
    def _from_asymmetric_unit(self, input_obj: AsymmetricUnit):
        BioAssemb._index+=1
        self.name = 'A.U.'
        self.transforms = TransformList.identity()

    @__init__.register(TransformList)
    def _from_transform_list(self, input_obj: TransformList):
        BioAssemb._index+=1
        self.name = f'Assembly{BioAssemb._index}'
        self.transforms = input_obj

    @__init__.register(List[PDBRecord])
    def _from_obj_list(self, input_obj: List[PDBRecord]):
        BioAssemb._index+=1
        self.name = f'Assembly{BioAssemb._index}'
        self.transforms = TransformList.from_pdb_records(input_obj)

    @classmethod
    def reset_index(cls):
        cls._index=1

    def activate(self, AU:AsymmetricUnit,CM:ChainIDManager):
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
            T.generate_chainIDmap(AU.segments.segnames, AU.segments.daughters, CM)

class BioAssembList(UserList):
    """
    A class for handling lists of BioAssemb objects.
    This class inherits from UserList and provides methods to manage
    collections of biological assemblies.
    """

    @singledispatchmethod
    def __init__(self, input_obj: Any):
        raise TypeError(f'Cannot initialize {type(self)} from object of type {type(input_obj)}')

    @__init__.register(Dict[str, PDBRecord])
    def _from_dict(self, input_obj: Dict[str, Any]):
        BioAssemb.reset_index()
        B=[]
        bareclabels=[x for x in input_obj if ('REMARK.350.BIOMOLECULE' in x and 'TRANSFORM' in x)]
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
            reclist=List([])
            savhdr=[]
            for t in trs:
                if f'REMARK.350.BIOMOLECULE{ba}.TRANSFORM{t}' in input_obj:
                    barec=input_obj[f'REMARK.350.BIOMOLECULE{ba}.TRANSFORM{t}']
                    if hasattr(barec,'header'):
                        savhdr=barec.header
                    else:
                        barec.header=savhdr
                    reclist.append(barec)
                    logger.debug(f'BA {ba} header {barec.header}')
                    logger.debug(barec.pstr())
            B.append(BioAssemb(reclist))
        logger.debug(f'There are {len(B)} biological assemblies')
        super().__init__(B)

    @__init__.register(DataContainer)
    def _from_data_container(self, input_obj: DataContainer):
        BioAssemb.reset_index()
        B=[]
        Assemblies=input_obj.getObj('pdbx_struct_assembly')
        gen=input_obj.getObj('pdbx_struct_assembly_gen')
        oper=input_obj.getObj('pdbx_struct_oper_list')
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
                    T=Transform.from_matvec(m,v,this_asyms,idx)
                    transforms.append(T)
                    idx+=1
        logger.debug(f'parsed {len(transforms)} transforms for ba {ba_idx}')
        BA=BioAssemb(transforms)
        B.append(BA)
        logger.debug(f'There '+'is' if len(B)==1 else 'are'+f' {len(B)} biological assembl'+'y' if len(B)==1 else 'ies')
        super().__init__(B)
