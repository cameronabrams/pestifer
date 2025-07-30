# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for processing biological assemblies
"""

import numpy as np
import re
# from functools import singledispatchmethod
from mmcif.api.PdbxContainers import DataContainer
from pidibble.pdbparse import get_symm_ops
from pidibble.pdbrecord import PDBRecord, PDBRecordList, PDBRecordDict
import logging
logger=logging.getLogger(__name__)
from typing import ClassVar
from collections import UserList

from .asymmetricunit import AsymmetricUnit
from .chainidmanager import ChainIDManager
from ..util.coord import build_tmat
from ..core.stringthings import plu

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
    def __init__(self, 
                 index: int = 0, 
                 tmat: np.ndarray | None = None, 
                 applies_chainIDs: list[str] | None = None, 
                 chainIDmap: dict[str, str] | None = None, 
                 segname_by_type_map: dict[str, str] | None = None):
        self.index = index
        self.tmat = tmat
        self.applies_chainIDs = applies_chainIDs
        self.chainIDmap = chainIDmap
        self.segname_by_type_map = segname_by_type_map

    @classmethod
    def from_pdb_record(cls, 
                        barec: PDBRecord, 
                        index: int = 0, 
                        applies_chainIDs: list[str] | None = None, 
                        segname_by_type_map: dict[str, str] | None = None):    
        RotMat, TransVec = get_symm_ops(barec)
        tmat = build_tmat(RotMat, TransVec)
        applies_chainIDs = barec.header if hasattr(barec, 'header') else []
        return cls(index=index, tmat=tmat, applies_chainIDs=applies_chainIDs, segname_by_type_map=segname_by_type_map)

    @classmethod
    def from_matvec(cls, 
                    RotMat: np.ndarray, TransVec: np.ndarray, 
                    index: int = 0, 
                    applies_chainIDs: list[str] | None = None, 
                    segname_by_type_map: dict[str, str] | None = None):
        tmat = build_tmat(RotMat, TransVec)
        return cls(index=index, tmat=tmat, applies_chainIDs=applies_chainIDs, segname_by_type_map=segname_by_type_map)

    @classmethod
    def identity(cls, 
                 index: int = 0,
                 applies_chainIDs: list[str] | None = None, 
                 segname_by_type_map: dict[str, str] | None = None):
        tmat = np.identity(4, dtype=float)
        applies_chainIDs = []
        return cls(index=index, tmat=tmat, applies_chainIDs=applies_chainIDs, segname_by_type_map=segname_by_type_map)

    def set_chainIDmap(self, chainIDmap: dict[str, str]):
        """
        Sets the chain ID mapping for the transformation.
        
        Parameters
        ----------
        chainIDmap : dict
            A dictionary mapping original chain IDs to new chain IDs.
        """
        self.chainIDmap = chainIDmap

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
            self.segname_by_type_map[segtype] = {}
        self.segname_by_type_map[segtype][chainID] = seglabel

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

    def generate_chainIDmap(self, auChainIDs: list[str], daughters: dict, CM: ChainIDManager):
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
            The ChainIDManager instance used to manage chain IDs.
        """
        applies_to=self.applies_chainIDs[:]
        for d,v in daughters.items():
            if d in self.applies_chainIDs:
                applies_to.extend(v)
        logger.debug(f'Transform applies to {applies_to}')
        if self.is_identity():
            logger.debug(f'Identity transform gets a thru map applied to {applies_to}')
            self.chainIDmap = CM.thru_map(auChainIDs, applies_to)
        else:
            logger.debug(f'Transform gets a new map applied to {applies_to}')
            self.chainIDmap = CM.generate_next_map(auChainIDs, applies_to)

class TransformList(UserList):
    """
    A class for handling lists of Transform objects.
    This class inherits from AncestorAwareObjList and provides methods to manage
    collections of Transform instances.
    """
    @classmethod
    def identity(cls):
        """
        Creates a TransformList containing a single identity transformation.
        """
        return cls([Transform('identity')])
            
    @classmethod
    def from_pdb_records(cls, pdb_records: PDBRecordList):
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

    def __init__(self, name: str | None = None, transforms: TransformList | None = None):
        """
        Initialize a BioAssemb instance.

        Parameters
        ----------
        name : str | None, optional
            The name of the biological assembly.
        transforms : TransformList | None, optional
            A list of Transform objects representing the transformations applied to the assembly.
        """
        self.name = name or f'Assembly{BioAssemb._index}'
        self.transforms = transforms or TransformList.identity()
        BioAssemb._index += 1

    @classmethod
    def from_dict(cls, input_obj: dict[str, str | TransformList | None]):
        if not ('name' in input_obj and 'transforms' in input_obj):
            raise ValueError(f'BioAssemb dictionary must contain "name" and "transforms" keys, got {input_obj}')
        if not isinstance(input_obj['transforms'], TransformList):
            raise TypeError(f'BioAssemb transforms must be a TransformList, got {type(input_obj["transforms"])}')
        cls._index += 1
        name = input_obj.get('name', f'Assembly{cls._index}')
        transforms = input_obj.get('transforms', TransformList.identity())
        return cls(name, transforms)

    @classmethod
    def from_asymmetric_unit(cls):
        """
        Initialize a BioAssemb instance from an AsymmetricUnit.
        
        Returns
        -------
        BioAssemb
            An instance of BioAssemb initialized with name A.U. and an identity transformation. 
        """
        cls._index += 1
        return cls(name='A.U.', transforms=TransformList.identity())

    @classmethod
    def from_transform_list(cls, TList: TransformList):
        """
        Initialize a BioAssemb instance from a TransformList.
        
        Parameters
        ----------
        input_obj : TransformList
            A TransformList containing the transformations for the biological assembly.
        
        Returns
        -------
        BioAssemb
            An instance of BioAssemb initialized with the provided TransformList.
        """
        cls._index += 1
        return cls(name=f'Assembly{cls._index}', transforms=TList)

    @classmethod
    def from_pdbrecordlist(cls, PRList: PDBRecordList):
        """
        Initialize a BioAssemb instance from a PDBRecordList.

        Parameters
        ----------
        PRList : PDBRecordList
            A PDBRecordList containing the transformations for the biological assembly.

        Returns
        -------
        BioAssemb
            An instance of BioAssemb initialized with the provided PDBRecordList.
        """
        cls._index += 1
        return cls(name=f'Assembly{cls._index}', transforms=TransformList.from_pdb_records(PRList))

    @classmethod
    def reset_index(cls):
        cls._index=1

    def activate(self, AU: AsymmetricUnit, CM: ChainIDManager):
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

    @classmethod
    def from_pdb_record_dict(cls, PRDict: PDBRecordDict):
        """
        Initialize a BioAssembList from a PDBRecordDict.
        
        Parameters
        ----------
        PRDict : PDBRecordDict
            A dictionary containing PDB records for biological assemblies.

        Returns
        -------
        BioAssembList
            An instance of BioAssembList initialized with the provided PDBRecordDict.
        """

        BioAssemb.reset_index()
        B=[]
        # search the PDBRecordDict for keys that match the pattern 'REMARK.350.BIOMOLECULE{n}.TRANSFORM{m}', and extract the assembly number n and transform number m
        # where n is the assembly number and m is the transform number
        ptn = r'REMARK\.350\.BIOMOLECULE(\d+)\.TRANSFORM(\d+)'
        savhdr = []
        records_of_ba: dict[int, list[PDBRecord]] = {}
        for key in PRDict.keys():
            match = re.match(ptn, key)
            if match:
                ba_number, transform_number = match.groups()
                ba_number  = int(ba_number)
                transform_number = int(transform_number)
                ba_record = PRDict[key]
                if hasattr(ba_record, 'header'):
                    savhdr = ba_record.header
                else:
                    ba_record.header = savhdr
                if not ba_number in records_of_ba:
                    records_of_ba[ba_number] = []
                records_of_ba[ba_number].append(ba_record)
        for ba_number, ba_recordlist in records_of_ba.items():
            logger.debug(f'BA {ba_number} has {len(ba_recordlist)} records')
            # Create a BioAssemb from the records
            B.append(BioAssemb.from_pdbrecordlist(ba_recordlist))
        logger.debug(f'There are {len(B)} biological assemblies')
        return cls(B)

    @classmethod
    def from_data_container(cls, dc: DataContainer):
        """
        Initialize a BioAssembList from a DataContainer.
        
        Parameters
        ----------
        dc : DataContainer
            A DataContainer containing the data for biological assemblies.

        Returns
        -------
        BioAssembList
            An instance of BioAssembList initialized with the provided DataContainer.
        """
        BioAssemb.reset_index()
        B = []
        Assemblies = dc.getObj('pdbx_struct_assembly')
        gen = dc.getObj('pdbx_struct_assembly_gen')
        oper = dc.getObj('pdbx_struct_oper_list')
        for ba_idx in range(len(Assemblies)):
            logger.debug(f'CIF: Establishing BA {ba_idx}')
            assemb_id = Assemblies.getValue('id', ba_idx)
            this_gen_idx_list = gen.selectIndices(assemb_id, 'assembly_id')
            logger.debug(f'BA {ba_idx} points to {len(this_gen_idx_list)} gen indexes')
            transforms = TransformList()
            for this_gen_idx in this_gen_idx_list:
                this_oper_list = gen.getValue('oper_expression', this_gen_idx).split(',')
                logger.debug(f'BA {ba_idx} gen {this_gen_idx} opers {this_oper_list}')
                this_asyms = gen.getValue('asym_id_list', this_gen_idx).split(',')
                logger.debug(f'asym ids: {this_asyms}')
                idx = 0
                # logger.debug(f'Expecting {len(this_opers)} transforms')
                for k, opere in enumerate(this_oper_list):
                    oper_idx = oper.selectIndices(opere, 'id')[0]
                    logger.debug(f'making transform from oper {oper_idx}')
                    m = np.identity(3)
                    v = np.zeros(3)
                    for i in range(3):
                        I = i + 1
                        vlabel = f'vector[{I}]'
                        v[i] = float(oper.getValue(vlabel, oper_idx))
                        for j in range(3):
                            J = j + 1
                            mlabel = f'matrix[{I}][{J}]'
                            m[i][j] = float(oper.getValue(mlabel, oper_idx))
                    T = Transform.from_matvec(m, v, this_asyms, idx)
                    transforms.append(T)
                    idx += 1
        logger.debug(f'parsed {len(transforms)} transforms for ba {ba_idx}')
        BA = BioAssemb(transforms)
        B.append(BA)
        logger.debug(f'There {plu(len(B), "is", "are")} {len(B)} biological assembl{plu(len(B), "y", "ies")}')
        return cls(B)
