# Author: Cameron F. Abrams, <cfa22@drexel.edu>

"""
A Transform class for handling homogeneous coordinate transformations of segments in a molecular structure.
This class represents a transformation that can be applied to segments in an asymmetric unit,
including rotation and translation.
"""
import logging
logger = logging.getLogger(__name__)
import numpy as np

from pidibble.pdbrecord import PDBRecord, PDBRecordList
from pidibble.pdbparse import get_symm_ops
from pydantic import Field
from typing import ClassVar

from .chainidmanager import ChainIDManager
from ..util.coord import build_tmat
from ..core.baseobj import BaseObj, BaseObjList

class Transform(BaseObj):
    """
    A class for handling transformations of segments in a molecular structure.
    This class represents a transformation that can be applied to segments
    in an asymmetric unit, including rotation and translation.

    This method used the :func:`pidibble.pdbparse.get_symm_ops` function to extract the rotation matrix and translation vector from a PDB record.
    It constructs a 4 x 4 transformation matrix from these components.
    """

    _required_fields = {'index', 'tmat', 'applies_chainIDs'}
    _optional_fields = {'chainIDmap', 'segname_by_type_map'}

    index: int = Field(..., description="Index of the transformation")
    tmat: np.ndarray = Field(..., description="4 x 4 transformation matrix")
    applies_chainIDs: list[str] = Field(..., description="List of chain IDs to which this transformation applies")
    chainIDmap: dict[str, str] | None = Field(None, description="Mapping of chain IDs for the transformation")
    segname_by_type_map: dict[str, str] | None = Field(None, description="Mapping of segment names by type for the transformation")

    _count: ClassVar[int] = 0

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for Transform instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        """

        if args and isinstance(args[0], PDBRecord):
            ba_record = args[0]
            RotMat, TransVec = get_symm_ops(ba_record)
            tmat = build_tmat(RotMat, TransVec)
            applies_chainIDs = ba_record.header if hasattr(ba_record, 'header') else []
            input_dict = dict(index=Transform._count, 
                              tmat=tmat, 
                              applies_chainIDs=applies_chainIDs)
            Transform._count += 1
            return input_dict
        elif args and isinstance(args[0], np.ndarray):
            rot = args[0]
            trans = args[1] if len(args) > 1 else np.zeros(3, dtype=float)
            if rot.shape != (3, 3):
                raise ValueError(f'Rotation matrix must be 3x3, got {rot.shape}')
            tmat = build_tmat(rot, trans)
            input_dict = dict(index=Transform._count, 
                              tmat=tmat, applies_chainIDs=[])
            Transform._count += 1
            return input_dict
        return super()._adapt(*args, **kwargs)
        
    @classmethod
    def identity(cls):
        """
        Creates an identity transformation.
        This method returns a Transform object with an identity transformation matrix and no chain ID mappings.
        """
        tmat = np.identity(4, dtype=float)
        input_dict = dict(index=cls._count, tmat=tmat, applies_chainIDs=[])
        cls._count += 1
        return cls(**input_dict)

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
        if not self.segname_by_type_map:
            self.segname_by_type_map = {}
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

class TransformList(BaseObjList[Transform]):
    """
    A class for handling lists of Transform objects.
    This class inherits from AncestorAwareObjList and provides methods to manage
    collections of Transform instances.
    """

    def describe(self):
        return f'<TransformList: {len(self)} items>'

    def __init__(self, initlist: list[Transform] | PDBRecordList = [], *args, **kwargs):
        if isinstance(initlist, PDBRecordList):
            # If the first argument is a PDBRecordList, convert it to a list of Transform objects
            initlist = [Transform(x) for x in initlist]
        super().__init__(initlist)

    @classmethod
    def identity(cls, length: int = 1):
        """
        Creates a TransformList containing a `length` identity transformations.
        """
        return cls([Transform.identity()] * length)

