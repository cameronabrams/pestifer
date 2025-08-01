# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A class for processing biological assemblies
"""

import logging
import numpy as np
import re
logger=logging.getLogger(__name__)

# from functools import singledispatchmethod
from mmcif.api.PdbxContainers import DataContainer
from pidibble.pdbrecord import PDBRecord, PDBRecordList, PDBRecordDict
from typing import ClassVar
from collections import UserList

from .asymmetricunit import AsymmetricUnit
from .chainidmanager import ChainIDManager
from ..core.stringthings import plu
from .transform import Transform, TransformList

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
