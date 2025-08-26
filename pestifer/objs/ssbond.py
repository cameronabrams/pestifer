# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Disufulfide bonds are covalent linkages between two CYS residues in a protein.
These are represented in PDB files as SSBOND records, and in mmCIF files
as struct_conn records.
"""
from __future__ import annotations

import logging

from mmcif.api.PdbxContainers import DataContainer
from pidibble.pdbrecord import PDBRecord, PDBRecordDict
from pydantic import Field
from typing import ClassVar, TYPE_CHECKING

from .mutation import MutationList
from .resid import ResID
from ..core.baseobj import BaseObj, BaseObjList
from ..psfutil.psfpatch import PSFDISUPatch
from ..util.cifutil import CIFdict

if TYPE_CHECKING:
    from ..molecule.residue import Residue, ResidueList

logger = logging.getLogger(__name__)

class SSBond(BaseObj):
    """
    A class for handling disulfide bonds between two CYS residues in a protein.
    """

    _required_fields = {'chainID1', 'resid1', 'chainID2', 'resid2'}
    """
    Required attributes for an SSBond object.
    These attributes must be provided when creating an SSBond object.
    
    - ``chainID1``: The chain ID of the first CYS residue.
    - ``resid1``: The residue ID of the first CYS residue.
    - ``chainID2``: The chain ID of the second CYS residue.
    - ``resid2``: The residue ID of the second CYS residue.
    """

    _optional_fields = {'serial_number', 'residue1', 'residue2', 'resname1', 'resname2', 'sym1', 'sym2', 'length', 'ptnr1_auth_asym_id', 'ptnr2_auth_asym_id', 'ptnr1_auth_seq_id', 'ptnr2_auth_seq_id'}
    """
    Optional attributes for an SSBond object.
    These attributes are not required but can be provided to enhance the SSBond object.
    
    - ``serial_number``: The serial number of the SSBond.
    - ``residue1``: The first CYS residue object associated with the SSBond.
    - ``residue2``: The second CYS residue object associated with the SSBond.
    - ``resname1``: The residue name of the first CYS residue.
    - ``resname2``: The residue name of the second CYS residue.
    - ``sym1``: The symmetry operator for the first CYS residue.
    - ``sym2``: The symmetry operator for the second CYS residue.
    - ``length``: The length of the SSBond.
    - ``ptnr1_auth_asym_id``: The author asym ID of the first CYS residue in mmCIF format.
    - ``ptnr2_auth_asym_id``: The author asym ID of the second CYS residue in mmCIF format.
    - ``ptnr1_auth_seq_id``: The author sequence ID of the first CYS residue in mmCIF format.
    - ``ptnr2_auth_seq_id``: The author sequence ID of the second CYS residue in mmCIF format.
    """

    chainID1: str = Field(..., description="Chain ID of the first CYS residue")
    resid1: ResID = Field(..., description="Residue ID of the first CYS residue")
    chainID2: str = Field(..., description="Chain ID of the second CYS residue")
    resid2: ResID = Field(..., description="Residue ID of the second CYS residue")

    serial_number: int = Field(0, description="Serial number of the SSBond")
    residue1: 'Residue' = Field(None, description="First CYS residue object associated with the SSBond")
    residue2: 'Residue' = Field(None, description="Second CYS residue object associated with the SSBond")
    resname1: str = Field('CYS', description="Residue name of the first CYS residue")
    resname2: str = Field('CYS', description="Residue name of the second CYS residue")
    sym1: str = Field('', description="Symmetry operator for the first CYS residue")
    sym2: str = Field('', description="Symmetry operator for the second CYS residue")
    length: float = Field(0.0, description="Length of the SSBond")
    ptnr1_auth_asym_id: str = Field('', description="Author asym ID of the first CYS residue in mmCIF format")
    ptnr2_auth_asym_id: str = Field('', description="Author asym ID of the second CYS residue in mmCIF format")
    ptnr1_auth_seq_id: str = Field('', description="Author sequence ID of the first CYS residue in mmCIF format")
    ptnr2_auth_seq_id: str = Field('', description="Author sequence ID of the second CYS residue in mmCIF format")

    _yaml_header: ClassVar[str] = 'ssbonds'
    """
    YAML header for SSBond objects.
    This header is used to identify SSBond objects in YAML files.
    """

    _objcat: ClassVar[str] = 'topol'
    """
    Category of the SSBond object.
    This categorization is used to group SSBond objects in the object manager.
    """

    _PDB_keyword: ClassVar[str] = 'SSBOND'
    """
    Keyword used to identify SSBond records in PDB files.
    """

    _CIF_CategoryName: ClassVar[str] = 'struct_conn'
    """
    Name of the CIF category that contains information for SSBond objects.
    """
    _CIF_CategoryElementTypes: ClassVar[dict[str, set]] = {'conn_type_id': {'disulf'}}
    """
    CIF category element types for SSBond objects; a CIF Category is effectively a list of dictionaries, and _CIF_CategoryElementTypes[keyname] is a set of values that are valid for that key, indicating the element is to be interpreted as a SSBond.
    """

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for SSBond instantiation.
        """
        if args and isinstance(args[0], str):
            # shortcode format: C_RRR-D_SSS
            # C, D chainIDs
            # RRR, SSS resids
            raw = args[0]
            s1, s2 = raw.split('-')
            resid1 = ResID(s1.split('_')[1])
            resid2 = ResID(s2.split('_')[1])
            return dict(chainID1=s1.split('_')[0], resid1=resid1, chainID2=s2.split('_')[0], resid2=resid2)
        elif args and isinstance(args[0], PDBRecord):
            pdbrecord = args[0]
            return dict(
                chainID1 = pdbrecord.residue1.chainID,
                resid1 = ResID(pdbrecord.residue1.seqNum, pdbrecord.residue1.iCode),
                chainID2 = pdbrecord.residue2.chainID,
                resid2 = ResID(pdbrecord.residue2.seqNum, pdbrecord.residue2.iCode),
                serial_number = pdbrecord.serNum,
                resname1 = 'CYS',
                resname2 = 'CYS',
                sym1 = pdbrecord.sym1,
                sym2 = pdbrecord.sym2,
                length = pdbrecord.length
            )
        elif args and isinstance(args[0], CIFdict):
            cd = args[0]
            return dict(
                chainID1 = cd['ptnr1_label_asym_id'],
                resid1 = ResID(cd['ptnr1_label_seq_id'], cd['pdbx_ptnr1_pdb_ins_code']),
                chainID2 = cd['ptnr2_label_asym_id'],
                resid2 = ResID(cd['ptnr2_label_seq_id'], cd['pdbx_ptnr2_pdb_ins_code']),
                serial_number = int(cd['id'].strip('disulf')),
                resname1 = 'CYS',
                resname2 = 'CYS',
                sym1 = cd['ptnr1_symmetry'],
                sym2 = cd['ptnr2_symmetry'],
                length = float(cd['pdbx_dist_value'])
            )
        elif args and isinstance(args[0], PSFDISUPatch):
            dp = args[0]
            return dict(
                chainID1 = dp.seg1,
                resid1 = dp.resid1,
                chainID2 = dp.seg2,
                resid2 = dp.resid2,
                # Default values for optional fields
                serial_number = 0,  # Default serial number
                resname1 = 'CYS',
                resname2 = 'CYS'
            )
        return super()._adapt(*args, **kwargs)

    def shortcode(self) -> str:
        """
        Convert the SSBond object to a shortcode string representation.
        
        Returns
        -------
        str
            A shortcode string representation of the SSBond object in the format C_RRR-D_SSS.
            Where C and D are chain IDs, and RRR and SSS are residue sequence numbers with insertion codes.
        """
        return f"{self.chainID1}_{self.resid1.resid}-{self.chainID2}_{self.resid2.resid}"

    def __str__(self):
        """
        Return a shortcode representation of the SSBond.
        The shortcode representation includes the chain IDs, residue sequence numbers, and insertion codes
        of both residues involved in the SSBond.
        If the residue attributes are set, it uses their chain IDs, sequence numbers, and insertion codes instead.
        """
        return self.shortcode()

    def pdb_line(self) -> str:
        """
        Write the SSBond as it would appear in a PDB file.

        Returns
        -------
        str
            The PDB line for the SSBOND record.
        """
        pdbline='SSBOND'+\
                '{:4d}'.format(self.serial_number)+\
                '{:>4s}'.format(self.resname1)+\
                '{:>2s}'.format(self.chainID1)+\
                '{:>5s}'.format(self.resid1.pdbresid)+\
                '{:>6s}'.format(self.resname2)+\
                '{:>2s}'.format(self.chainID2)+\
                '{:>5s}'.format(self.resid2.pdbresid)+\
                ' '*23+\
                '{:>6s}'.format(self.sym1)+\
                '{:>7s}'.format(self.sym2)+\
                '{:6.2f}'.format(self.length)
        return pdbline

class SSBondList(BaseObjList[SSBond]):
    """
    A class for handling a list of SSBond objects.
    This class inherits from BaseObjList and provides methods to manage
    a list of SSBond objects.
    """

    def describe(self):
        return f'SSBondList with {len(self)} SSBonds'

    @classmethod
    def from_pdb(cls, pdb: PDBRecordDict) -> SSBondList:
        """
        Create a SSBondList from a PDBRecordDict.
        """
        if SSBond._PDB_keyword not in pdb:
            return cls([])
        return cls([SSBond(x) for x in pdb[SSBond._PDB_keyword]])

    @classmethod
    def from_cif(cls, dc: DataContainer) -> SSBondList:
        """
        Create a SSBondList from a CIF DataContainer.

        Parameters
        ----------
        dc : DataContainer
            A CIF DataContainer containing the necessary fields to create SSBond objects.

        Returns
        -------
        SSBondList
            An instance of SSBondList created from the CIF DataContainer.
        """
        L = []
        cif_category = dc.getObj(SSBond._CIF_CategoryName)
        if cif_category is None:
            return cls([])
        for i in range(len(cif_category)):
            for key, valset in SSBond._CIF_CategoryElementTypes.items():
                objTypeid = cif_category.getValue(key, i)
                if objTypeid in valset:
                    this_link = SSBond(CIFdict(cif_category, i))
                    L.append(this_link)
        return cls(L)

    def assign_residues(self, Residues: 'ResidueList'):
        """
        Assigns a list of Residue objects to the SSBond residues.
        
        Parameters
        ----------
        Residues : list of Residue
            A list of Residue objects to be assigned to the SSBond residues.
        
        Returns
        -------
        SSBondList
            A new SSBondList object containing the residues that were not assigned because no applicable residues were found.
        """
        logger.debug(f'SSBonds: Assigning residues from list of {len(Residues)} residues')
        ignored_by_ptnr1=self.assign_objs_to_attr('residue1', Residues, chainID='chainID1', resid='resid1')
        ignored_by_ptnr2=self.assign_objs_to_attr('residue2', Residues, chainID='chainID2', resid='resid2')
        return self.__class__(ignored_by_ptnr1 + ignored_by_ptnr2)

    # def prune_mutations(self, Mutations: MutationList) -> 'SSBondList':
    #     """
    #     Prunes the SSBondList by removing SSBonds that are associated with mutations. This method iterates over each Mutation in the provided list and checks if there are
    #     any SSBonds that match the chain ID, residue sequence number, and insertion code of
    #     the Mutation. If a match is found, the SSBond is removed from the list 
    #     and added to a new SSBondList object. The method returns this new SSBondList object.
        
    #     Parameters
    #     ----------
    #     Mutations : MutationList
    #         A list of Mutation objects that contain information about the mutations to be pruned.
        
    #     Returns
    #     -------
    #     SSBondList
    #         A new SSBondList object containing the SSBonds that were pruned.

    #     """
    #     pruned = self.__class__([])
    #     for m in Mutations.data:
    #         left = self.get(chainID1=m.chainID, resid1=m.resid)
    #         if left:
    #             pruned.take_item(left, self)
    #         right = self.get(chainID2=m.chainID, resid2=m.resid)
    #         if right:
    #             pruned.take_item(right, self)
    #     return pruned

    # def map_attr(self, mapped_attr: str, key_attr: str, map: dict):
    #     """
    #     Maps an attribute of the SSBondList to a new value using a mapping dictionary.  A dictionary that maps the values of the key_attr to new values for the mapped_attr.
    #     This method iterates over each SSBond in the list and applies the mapping to the
    #     specified attribute. It logs the mapping operation and the number of SSBonds being processed.
        
    #     Parameters
    #     ----------
    #     mapped_attr : str
    #         The attribute of the SSBondList to be mapped.
    #     key_attr : str
    #         The attribute of the SSBond that will be used as the key for the mapping.
    #     map : dict
    #     """
    #     logger.debug(f'Mapping {mapped_attr} {key_attr} in {len(self)} SSBonds using map {map}')
    #     super().map_attr(mapped_attr, key_attr, map)
