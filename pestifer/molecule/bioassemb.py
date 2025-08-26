# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for processing biological assemblies
"""
from __future__ import annotations

import logging
import re

import numpy as np

from mmcif.api.PdbxContainers import DataContainer
from pidibble.pdbrecord import PDBRecord, PDBRecordList, PDBRecordDict
from typing import ClassVar, TYPE_CHECKING
from collections import UserList

from .asymmetricunit import AsymmetricUnit
from .chainidmanager import ChainIDManager
if TYPE_CHECKING:
    from .molecule import Molecule
from .transform import Transform, TransformList

from ..util.stringthings import plu

logger = logging.getLogger(__name__)

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

    def __init__(self, *args, **kwargs):
        """
        Initialize a BioAssemb instance.

        Parameters
        ----------
        name : str | None, optional
            The name of the biological assembly.
        transforms : TransformList | None, optional
            A list of Transform objects representing the transformations applied to the assembly.
        """
        apparent_index = kwargs.get('index', 0)
        if isinstance(apparent_index, str):
            if not apparent_index.isdigit():
                raise ValueError(f'Index must be an integer or a string representing an integer, got {apparent_index}')
            apparent_index = int(apparent_index)
        self.index = apparent_index
        self.name = kwargs.get('name', f'Assembly{self.index}')
        logger.debug(f'Initializing BioAssemb {self.name} with index {self.index}')
        if len(args) == 0:
            if len(kwargs) == 0:
                self.transforms = TransformList.identity(1)
            else:
                if 'transforms' in kwargs:
                    self.transforms = kwargs.get('transforms', TransformList.identity(1))
                elif 'pdbrecordlist' in kwargs:
                    self.transforms = TransformList(kwargs['pdbrecordlist'])
                elif 'asymmetric_unit' in kwargs:
                    self.transforms = TransformList.identity(1)
                else:
                    raise ValueError(f'BioAssemb requires either a name or transforms, got {kwargs}')           
        elif len(args) == 1:
            if isinstance(args[0], TransformList):
                self.transforms = args[0]
            elif isinstance(args[0], PDBRecordList):
                self.transforms = TransformList(args[0])
            elif isinstance(args[0], AsymmetricUnit):
                self.transforms = TransformList.identity(1)
            elif isinstance(args[0], dict):
                self.transforms = args[0].get('transforms', TransformList.identity(1))
            else:
                logger.warning(f'I do not know how to initialize BioAssemb with {args[0]} (type{type(args[0])})')
        elif len(args) == 2:
            self.name = args[0]
            self.transforms = args[1]
        else:
            logger.warning(f'I really do not know how to initialize BioAssemb with {args[0]} (type{type(args[0])})')
        self.parent_molecule = None

    def set_parent_molecule(self, mol: 'Molecule'):
        """
        Set the parent molecule for this BioAssemb instance.
        
        Parameters
        ----------
        mol : Molecule
            The parent molecule to set.
        """
        self.parent_molecule = mol

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
        for T in self.transforms.data:
            T.generate_chainIDmap(AU.segments.segnames, AU.segments.daughters, CM)

class BioAssembList(UserList[BioAssemb]):
    """
    A class for handling lists of BioAssemb objects.
    This class inherits from UserList and provides methods to manage
    collections of biological assemblies.
    """

    def __init__(self, initial_data=None):
        """
        Initialize a BioAssembList instance.
        
        Parameters
        ----------
        initial_data : list[BioAssemb] | None, optional
            Initial data for the BioAssembList. If None, an empty list is created.
        """
        if isinstance(initial_data, BioAssemb):
            initial_data = [initial_data]
        elif initial_data is None:
            initial_data = []
        elif isinstance(initial_data, PDBRecordDict):
            initial_data = BioAssembList.from_pdb_record_dict(initial_data)
        elif isinstance(initial_data, DataContainer):
            initial_data = BioAssembList.from_data_container(initial_data)
        super().__init__(initial_data or [])
        self.parent_molecule = None

    def get(self, index: int) -> BioAssemb | None:
        """
        Get a BioAssemb instance by index.

        Parameters
        ----------
        index : int
            The index of the BioAssemb to retrieve.

        Returns
        -------
        BioAssemb
            The BioAssemb instance at the specified index.
        """
        for ba in self.data:
            if ba.index == index:
                return ba
        return None

    def set_parent_molecule(self, mol: 'Molecule'):
        """
        Set the parent molecule for this BioAssembList instance.
        
        Parameters
        ----------
        mol : object
            The parent molecule to set.
        """
        self.parent_molecule = mol
        for ba in self.data:
            ba.set_parent_molecule(mol)

    @staticmethod
    def from_pdb_record_dict(PRDict: PDBRecordDict):
        """
        Initialize a BioAssembList from a PDBRecordDict.
        
        Parameters
        ----------
        PRDict : PDBRecordDict
            A dictionary containing PDB records for biological assemblies.

        Returns
        -------
        list
            An instance of BioAssembList initialized with the provided PDBRecordDict.
        """

        B=[]
        # search the PDBRecordDict for keys that match the pattern 'REMARK.350.BIOMOLECULE{n}.TRANSFORM{m}', and extract the assembly number n and transform number m
        # where n is the assembly number and m is the transform number
        ptn = r'REMARK\.350\.BIOMOLECULE(\d+)\.TRANSFORM(\d+)'
        savhdr = []
        records_of_ba: dict[int, PDBRecordList] = {}
        for key in PRDict.keys():
            match = re.match(ptn, key)
            if match:
                ba_number, transform_number = match.groups()
                ba_number = int(ba_number)
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
            B.append(BioAssemb(PDBRecordList(ba_recordlist), index=ba_number))
        logger.debug(f'There {plu(len(B), "is", "are")} {len(B)} biological assembl{plu(len(B), "y", "ies")}')
        return B

    @staticmethod
    def from_data_container(dc: DataContainer):
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
                    T = Transform(m, v, this_asyms, idx)
                    transforms.append(T)
                    idx += 1
            logger.debug(f'parsed {len(transforms)} transforms for assemb_id {assemb_id}')
            BA = BioAssemb(transforms, index=assemb_id)
            B.append(BA)
        logger.debug(f'There {plu(len(B), "is", "are")} {len(B)} biological assembl{plu(len(B), "y", "ies")}')
        return B
