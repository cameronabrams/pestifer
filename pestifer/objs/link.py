# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A Link is a covalent bond between two residues in a protein structure.
"""
from __future__ import annotations

import logging

import numpy as np

from mmcif.api.PdbxContainers import DataContainer
from pidibble.pdbrecord import PDBRecord, PDBRecordDict
from pydantic import Field
from typing import ClassVar, TYPE_CHECKING

from .mutation import MutationList
from .resid import ResID

from ..core.baseobj import BaseObj, BaseObjList
from ..util.coord import measure_dihedral

if TYPE_CHECKING:
    from ..molecule.atom import Atom
    from ..molecule.residue import Residue, ResidueList
    from ..molecule.segment import Segment, SegmentList

from ..psfutil.psfpatch import PSFLinkPatch

from ..util.cifutil import CIFdict

logger = logging.getLogger(__name__)

class Link(BaseObj):
    """
    A class for handling covalent bonds between residues where at least one residue is non-protein
    """
    _required_fields = {'chainID1', 'resid1', 'name1',
                        'chainID2', 'resid2', 'name2'}

    _optional_fields = {'altloc1', 'altloc2', 'resname1', 'resname2', 'sym1', 'sym2', 'link_distance', 'segname1', 'segname2', 'residue1', 'residue2', 'atom1', 'atom2', 'empty', 'segtype1', 'segtype2', 'ptnr1_label_asym_id', 'ptnr2_label_asym_id', 'ptnr1_label_seq_id', 'ptnr2_label_seq_id', 'ptnr1_label_comp_id', 'ptnr2_label_comp_id', 'ptnr1_auth_asym_id', 'ptnr2_auth_asym_id', 'ptnr1_auth_seq_id', 'ptnr2_auth_seq_id', 'ptnr1_auth_comp_id', 'ptnr2_auth_comp_id','patchname','patchhead'}

    _attr_choices = {'patchhead': {1, 2}}

    chainID1: str = Field(..., description="Chain ID of the first residue in the link")
    resid1: ResID = Field(..., description="Residue ID of the first residue in the link")
    name1: str = Field(..., description="Name of the first atom in the link")
    chainID2: str = Field(..., description="Chain ID of the second residue in the link")
    resid2: ResID = Field(..., description="Residue ID of the second residue in the link")
    name2: str = Field(..., description="Name of the second atom in the link")
    """
    Required attributes for a Link object.
    These attributes must be provided when creating a Link object.

    - ``chainID1``: The chain ID of the first residue in the link.
    - ``resid1``: The residue ID of the first residue in the link.
    - ``name1``: The name of the first atom in the link.
    - ``chainID2``: The chain ID of the second residue in the link.
    - ``resid2``: The residue ID of the second residue in the link.
    - ``name2``: The name of the second atom in the link.
    """
    altloc1: str | None = Field(None, description="Alternate location identifier for the first atom")
    altloc2: str | None = Field(None, description="Alternate location identifier for the second atom")
    resname1: str | None = Field(None, description="Residue name of the first residue in the link")
    resname2: str | None = Field(None, description="Residue name of the second residue in the link")
    sym1: str | None = Field(None, description="Symmetry operator for the first residue")
    sym2: str | None = Field(None, description="Symmetry operator for the second residue")
    link_distance: float | None = Field(None, description="Distance between the two atoms in the link")
    segname1: str | None = Field(None, description="Segment name of the first residue")
    segname2: str | None = Field(None, description="Segment name of the second residue")
    residue1: "Residue" = Field(None, description="First residue object in the link")
    residue2: "Residue" = Field(None, description="Second residue object in the link")
    atom1: "Atom" = Field(None, description="First atom object in the link")
    atom2: "Atom" = Field(None, description="Second atom object in the link")
    empty: bool | None = Field(None, description="Indicates if the link is empty")
    segtype1: str | None = Field(None, description="Segment type of the first residue")
    segtype2: str | None = Field(None, description="Segment type of the second residue")
    ptnr1_label_asym_id: str | None = Field(None, description="Asym ID of the first partner in the link (mmCIF)")
    ptnr2_label_asym_id: str | None = Field(None, description="Asym ID of the second partner in the link (mmCIF)")
    ptnr1_label_seq_id: str | None = Field(None, description="Sequence ID of the first partner in the link (mmCIF)")
    ptnr2_label_seq_id: str | None = Field(None, description="Sequence ID of the second partner in the link (mmCIF)")
    ptnr1_label_comp_id: str | None = Field(None, description="Component ID of the first partner in the link (mmCIF)")
    ptnr2_label_comp_id: str | None = Field(None, description="Component ID of the second partner in the link (mmCIF)")
    ptnr1_auth_asym_id: str | None = Field(None, description="Author asym ID of the first partner in the link (mmCIF)")
    ptnr2_auth_asym_id: str | None = Field(None, description="Author asym ID of the second partner in the link (mmCIF)")
    ptnr1_auth_seq_id: str | None = Field(None, description="Author sequence ID of the first partner in the link (mmCIF)")
    ptnr2_auth_seq_id: str | None = Field(None, description="Author sequence ID of the second partner in the link (mmCIF)")
    ptnr1_auth_comp_id: str | None = Field(None, description="Author component ID of the first partner in the link (mmCIF)")
    ptnr2_auth_comp_id: str | None = Field(None, description="Author component ID of the second partner in the link (mmCIF)")
    patchname: str | None = Field(None, description="Name of the patch applied to the link")
    patchhead: int | None = Field(None, description="1 = residue1 is the first residue in the link, 2 = residue2 is the first residue in the link")

    """
    Optional attributes for a Link object.
    These attributes can be provided to modify the behavior of the link.
    - ``residue1``: The first residue object in the link.
    - ``residue2``: The second residue object in the link.
    - ``atom1``: The first atom object in the link.
    - ``atom2``: The second atom object in the link.
    - ``empty``: A boolean indicating if the link is empty.
    - ``segtype1``: The segment type of the first residue.
    - ``segtype2``: The segment type of the second residue.
    - ``ptnr1_label_asym_id``: The asym ID of the first partner in the link (mmCIF).
    - ``ptnr2_label_asym_id``: The asym ID of the second partner in the link (mmCIF).
    - ``ptnr1_label_seq_id``: The sequence ID of the first partner in the link (mmCIF).
    - ``ptnr2_label_seq_id``: The sequence ID of the second partner in the link (mmCIF).
    - ``ptnr1_label_comp_id``: The component ID of the first partner in the link (mmCIF).
    - ``ptnr2_label_comp_id``: The component ID of the second partner in the link (mmCIF).
    - ``ptnr1_auth_asym_id``: The author asym ID of the first partner in the link (mmCIF).
    - ``ptnr2_auth_asym_id``: The author asym ID of the second partner in the link (mmCIF).
    - ``ptnr1_auth_seq_id``: The author sequence ID of the first partner in the link (mmCIF).
    - ``ptnr2_auth_seq_id``: The author sequence ID of the second partner in the link (mmCIF).
    - ``ptnr1_auth_comp_id``: The author component ID of the first partner in the link (mmCIF).
    - ``ptnr2_auth_comp_id``: The author component ID of the second partner in the link (mmCIF).
    - ``patchname``: The name of the patch applied to the link.
    - ``patchhead``: 1 = residue1 is the first residue in the link, 2 = residue2 is the first residue in the link.
    """
    
    _yaml_header: ClassVar[str] = 'links'
    """
    YAML header for Link objects.
    This header is used to identify Link objects in YAML files.
    """

    _PDB_keyword: ClassVar[str] = 'LINK'
    """
    PDB keyword for Link objects.
    """

    _CIF_CategoryName: ClassVar[str] = 'struct_conn'
    """
    Name of the CIF category that contains information for Link objects.
    """
    _CIF_CategoryElementTypes: ClassVar[dict[str, set]] = {'conn_type_id': {'covale', 'metalc'}}
    """
    CIF category element types for Link objects; a CIF Category is effectively a list of dictionaries, and _CIF_CategoryElementTypes[keyname] is a set of values that are valid for that key, indicating the element is to be interpreted as a Link.
    """

    _objcat: ClassVar[str] = 'topol'
    """
    Category of the Link object.
    This categorization is used to group Link objects in the object manager.
    """

    _patch_atomnames: ClassVar[dict[str, list[str]]] = {
        'NGLA':['ND2','C1'],
        'NGLB':['ND2','C1'],
        'SGPA':['OG','C1'],
        'SGPB':['OG','C1'],
        '11aa':['O1','C1'],
        '11ab':['O1','C1'],
        '11bb':['O1','C1'],
        '12aa':['O2','C1'],
        '12ab':['O2','C1'],
        '12ba':['O2','C1'],
        '12bb':['O2','C1'],
        '13aa':['O3','C1'],
        '13ab':['O3','C1'],
        '13ba':['O3','C1'],
        '13bb':['O3','C1'],
        '14aa':['O4','C1'],
        '14ab':['O4','C1'],
        '14ba':['O4','C1'],
        '14bb':['O4','C1'],
        '16AT':['O6','C1'],
        '16BT':['O6','C1'],
        'SA26AT':['O6','C2'],
        'SA28AA':['O8','C2'],
        'SA29AT':['O9','C2'],
        'ZNHE':['ZN','NE2'],
        'ZNHD':['ZN','ND2']
    }
    """
    A dictionary mapping atom names to their corresponding CHARMM36 patch names.
    """

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for Link instantiation.
        This method is used to convert various input types into a dictionary of parameters.
        """
        if args and isinstance(args[0], str):
            input_dict = Link._from_shortcode(args[0])
            return input_dict
        elif args and isinstance(args[0], PDBRecord):
            return Link._from_pdbrecord(args[0])
        elif args and isinstance(args[0], CIFdict):
            return Link._from_cifdict(args[0])
        elif args and isinstance(args[0], PSFLinkPatch):
            return Link._from_psflinkpatch(args[0])
        return super()._adapt(*args, **kwargs)
    
    @staticmethod
    def _from_pdbrecord(pdbrecord: PDBRecord) -> dict:
        return {
            'name1': pdbrecord.name1,
            'resname1': pdbrecord.residue1.resName,
            'chainID1': pdbrecord.residue1.chainID,
            'resid1': ResID(pdbrecord.residue1.seqNum, pdbrecord.residue1.iCode),
            'name2': pdbrecord.name2,
            'resname2': pdbrecord.residue2.resName,
            'chainID2': pdbrecord.residue2.chainID,
            'resid2': ResID(pdbrecord.residue2.seqNum, pdbrecord.residue2.iCode),
            'altloc2': pdbrecord.altLoc2,
            'altloc1': pdbrecord.altLoc1,
            'sym1': pdbrecord.sym1,
            'sym2': pdbrecord.sym2,
            'link_distance': pdbrecord.length,
            'segname1': pdbrecord.residue1.chainID,
            'segname2': pdbrecord.residue2.chainID,
            'empty': False
        }

    @staticmethod
    def _from_cifdict(cd: CIFdict) -> dict:
        resseqnum1 = int(cd['ptnr1_label_seq_id']) if cd['ptnr1_label_seq_id'] != '.' else int(cd['ptnr1_auth_seq_id'])
        resseqnum2 = int(cd['ptnr2_label_seq_id']) if cd['ptnr2_label_seq_id'] != '.' else int(cd['ptnr2_auth_seq_id'])
        insertion1 = cd['pdbx_ptnr1_pdb_ins_code'] if cd['pdbx_ptnr1_pdb_ins_code'] != '.' else ''
        insertion2 = cd['pdbx_ptnr2_pdb_ins_code'] if cd['pdbx_ptnr2_pdb_ins_code'] != '.' else ''
        resid1 = ResID(resseqnum1, insertion1)
        resid2 = ResID(resseqnum2, insertion2)
        return {
                'name1': cd['ptnr1_label_atom_id'],
                'altloc1': cd['pdbx_ptnr1_label_alt_id'],
                'resname1': cd['ptnr1_label_comp_id'],
                'chainID1': cd['ptnr1_label_asym_id'],
                'resid1': resid1,
                'name2': cd['ptnr2_label_atom_id'],
                'altloc2': cd['pdbx_ptnr2_label_alt_id'],
                'resname2': cd['ptnr2_label_comp_id'],
                'chainID2': cd['ptnr2_label_asym_id'],
                'resid2': resid2,
                'sym1': cd.get('ptnr1_symmetry', ''),
                'sym2': cd.get('ptnr2_symmetry', ''),
                'link_distance': float(cd.get('pdbx_dist_value', 0.0)),
                'segname1': cd['ptnr1_label_asym_id'],
                'segname2': cd['ptnr2_label_asym_id'],
                'empty': False
            }

    @staticmethod
    def _from_psflinkpatch(L: PSFLinkPatch) -> dict:
        idict = {
            'chainID1': L.seg1,
            'resid1': L.resid1,
            'chainID2': L.seg2,
            'resid2': L.resid2,
            'patchname': L.patchname if hasattr(L, 'patchname') else '',
            'patchhead': 1,  # default order
            'name1': Link._patch_atomnames[L.patchname][0],
            'name2': Link._patch_atomnames[L.patchname][1],
            'empty': False,
            }
        return idict

    @staticmethod
    def _from_shortcode(raw: str) -> dict:
        """
        Create a Link.Adapter instance from a string representation of a link.
        The string should be in the format 'C1_R1_A1-C2_R2_A2', where:
        - C1 is the chain ID of the first residue
        - R1 is the residue ID of the first residue
        - A1 is the atom name of the first residue that is part of the link
        - C2 is the chain ID of the second residue
        - R2 is the residue ID of the second residue
        - A2 is the atom name of the second residue that is linked to A1 of R1 of C1
        """
        I, J = raw.split('-')
        s1, ri1, a1 = I.split('_')
        resid1 = ResID(ri1)
        s2, ri2, a2 = J.split('_')
        resid2 = ResID(ri2)
        input_dict = {
            'chainID1': s1,
            'resid1': resid1,
            'name1': a1,
            'chainID2': s2,
            'resid2': resid2,
            'name2': a2,
            'empty': False
        }
        return input_dict

    def shortcode(self) -> str:
        """
        Returns a string representation of the link in the format 'C1_R1_A1-C2_R2_A2'.
        """
        return f"{self.chainID1}_{self.resid1.resid}_{self.name1}-{self.chainID2}_{self.resid2.resid}_{self.name2}"

    def set_patchname(self, force=False):
        """
        Set the charmmff patch name for this link.
        This method assigns a patch name based on the residues involved in the link.
        If the patch name is already set and ``force`` is False, it will not change the patch name.
        This method does not return any value. It modifies the ``patchname`` attribute of the Link object.
        It checks the residues involved in the link and assigns a patch name based on predefined mappings.
        
        Parameters
        ----------
        force : bool, optional
            If True, forces the patch name to be set even if it is already assigned.
            Default is False.
        """
        if hasattr(self,'patchname') and (self.patchname and len(self.patchname)>0) and not force:
            logger.debug(f'Patchname for {str(self)} already set to {self.patchname}')
            return
        self.patchname=''
        self.patchhead=1
        logger.debug(f'patch assignment for link {str(self)}')
        if not self.residue1 and not self.residue2:
            logger.debug(f'missing residue')
            logger.debug(f'1 {self.residue1}')
            logger.debug(f'2 {self.residue2}')
            return
        logger.debug(f'resname1 {self.resname1} segtype2 {self.segtype2}')
        my_res12 = [self.residue1, self.residue2]
        if self.resname1 == 'ASN' and self.segtype2 == 'glycan':
            # N-linked glycosylation site (toppar_all36_carb_glycopeptide)
            ICmap=[
                {'ICatomnames': ['1CG','1ND2','2C1','2O5'],
                 'mapping': {'NGLA':168.99,'NGLB':-70.91}}
            ]
            self.patchname = ic_reference_closest(my_res12, ICmap)
        elif self.resname1 == 'SER' and self.segtype2 == 'glycan':
            # O-linked to serine (toppar_all36_carb_glycopeptide)
            ICmap=[
                {'ICatomnames':['1CB','1OG','2C1','2O5'],
                 'mapping':{'SGPA':45.37,'SGPB':19.87}}
            ]
            self.patchname = ic_reference_closest(my_res12, ICmap)
        elif self.resname1 == 'THR' and self.segtype2 == 'glycan':
            # O-linked to threonine (toppar_all36_carb_glycopeptide)
            ICmap=[
                {'ICatomnames':['1CB','1OG1','2C1','2O5'],
                 'mapping':{'SGPA':69.9,'SGPB':33.16}}
            ]
            self.patchname = ic_reference_closest(my_res12, ICmap)
        elif self.name2 == 'C1' and self.segtype2 == 'glycan' and self.segtype1 == 'glycan':
            # all taken from 1xyz pres in top_all36_carb.rtf
            # including PHI angles of ICs with atoms in both residues
            if self.name1 == 'O1':  # 1->1 link
                ICmap=[
                    {'ICatomnames':'1O5  1C1  1O1  2C1'.split(),
                     'mapping':{'11aa':103.46,'11ab':121.75,'11bb':-56.58}
                    },
                    {'ICatomnames':'1C1  1O1  2C1  2O5'.split(),
                     'mapping':{'11aa':103.54,'11ab':51.80,'11bb':-79.64}
                    },
                    {'atomnames':'1O1  2C1  2O5  2C5'.split(),
                     'mapping':{'11aa':64.56,'11ab':167.51,'11bb':172.18}}
                ]
                self.patchname = ic_reference_closest(my_res12, ICmap)
            elif self.name1 == 'O2':  # 1->2 link
                ICmap=[
                    {'ICatomnames':'1C1  1C2  1O2  2C1'.split(),
                     'mapping':{'12aa':-132.81,'12ab':115.32,'12ba':-133.78,'12bb':117.14}
                    },
                    {'ICatomnames':'1C2  1O2  2C1  2O5'.split(),
                     'mapping':{'12aa':47.16,'12ab':86.93,'12ba':168.07,'12bb':-168.07}
                    }
                ]
                self.patchname = ic_reference_closest(my_res12, ICmap)
            elif self.name1 == 'O3':  # 1->3 link
                ICmap=[
                    {'ICatomnames':'1C2  1C3  1O3  2C1'.split(),
                     'mapping':{'13aa':113.19,'13ab':-141.32,'13ba':-131.68,'13bb':-141.32}
                    },
                    {'ICatomnames':'1C3  1O3  2C1  2O5'.split(),
                     'mapping':{'13aa':65.46,'13ab':65.46,'13ba':-100.16,'13bb':-130.16}
                    }
                ]
                self.patchname = ic_reference_closest(my_res12, ICmap)
            elif self.name1 == 'O4':  # 1->4 link
                ICmap=[
                    {'ICatomnames':'1C3  1C4  1O4  2C1'.split(),
                     'mapping':{'14aa':-86.29,'14ab':72.71,'14ba':-86.3,'14bb':81.86}},
                    {'ICatomnames':'1C4  1O4  2C1  2O5'.split(),
                     'mapping':{'14aa':133.57,'14ab':48.64,'14ba':-130.97,'14bb':-130.97}}
                ]
                self.patchname = ic_reference_closest(my_res12, ICmap)
            elif self.name1 == 'O6':  # 1->6 link
                ICmap=[
                    {'ICatomnames':'1C6  1O6  2C1  2O5'.split(),
                     'mapping':{'16AT':71.24,'16BT':-63.49}}
                ]
                self.patchname = ic_reference_closest(my_res12, ICmap)
        elif self.name2 == 'C2' and self.segtype2 == 'glycan' and self.segtype1 == 'glycan':
            if self.name1 == 'O6':
                self.patchname = 'SA26AT'
            elif self.name1 == 'O8':
                self.patchname = 'SA28AA'
            elif self.name1 == 'O9':
                self.patchname = 'SA29AT'
        elif self.name1 == 'O6' and self.name2 == 'C2':
            self.patchname = 'SA26AT'
        elif 'ZN' in self.resname1 and 'HIS' in self.resname2:
                if self.name2 == 'NE2':
                    self.patchname = 'ZNHE'
                else:
                    self.patchname = 'ZNHD'
                self.patchhead = 2
        elif 'HIS' in self.resname1 and 'ZN' in self.resname2:
                if self.name1 == 'NE2':
                    self.patchname = 'ZNHE'
                else:
                    self.patchname = 'ZNHD'
        elif 'HIS' in self.resname1 and 'HEM' in self.resname2:
                if self.name1 == 'NE2':
                    self.patchname = 'PHEM'
                else:
                    self.patchname = 'UNFOUND'
        elif 'HEM' in self.resname1 and 'HIS' in self.resname2:
                if self.name2 == 'NE2':
                    self.patchname = 'PHEM'
                else:
                    self.patchname = 'UNFOUND'
                self.patchhead = 2
        else:
            logger.warning(f'Could not identify patch for link {self.resname1}-{self.resname2}')
            self.patchname = 'UNFOUND'

    def update_residue(self, idx, **fields):
        """
        Updates the chainID of the residue in the link based on the index provided
        
        Parameters
        ----------
        idx: int
            1 for residue1, 2 for residue2
        fields: dict
            a dictionary of fields to update; currently only 'chainID' is supported
        """
        if idx==1:
            if 'chainID' in fields:
                new_chainID=fields['chainID']
                if self.chainID1!=new_chainID:
                    logger.debug(f'updating link residue1 chainID from {self.chainID1} to {new_chainID}')
                    self.chainID1=new_chainID
                    self.residue1.set(chainID=new_chainID)
        elif idx==2:
            if 'chainID' in fields:
                new_chainID=fields['chainID']
                if self.chainID2!=new_chainID:
                    logger.debug(f'updating link residue2 chainID from {self.chainID2} to {new_chainID}')
                    self.chainID2=new_chainID
                    self.residue2.set(chainID=new_chainID)

    def __str__(self):
        # if this a nascent link object that hasn't been processed by assigning residue objects to it,
        # just regurgitate the chainID and residue ID information that was used to create it.
        if not hasattr(self,'residue1') or not self.residue1 or not hasattr(self,'residue2') or not self.residue2:
            return f'{self.chainID1}_{self.resname1}{self.resid1.resid}-{self.chainID2}_{self.resname2}{self.resid2.resid}'
        # otherwise, access the string representation of the residues
        return f'{str(self.residue1)}-{str(self.residue2)}'

class LinkList(BaseObjList[Link]):
    """
    A class for handling lists of Links
    """

    @classmethod
    def from_pdb(cls, pdb: PDBRecordDict) -> 'LinkList':
        """
        Create a LinkList from a PDBRecordDict.
        """
        if Link._PDB_keyword not in pdb:
            return cls([])
        return cls([Link(x) for x in pdb[Link._PDB_keyword]])

    @classmethod
    def from_cif(cls, dc: DataContainer) -> 'LinkList':
        """
        Create a LinkList from a CIF DataContainer.
        
        Parameters
        ----------
        dc : DataContainer
            A CIF DataContainer containing the necessary fields to create Link objects.
        
        Returns
        -------
        LinkList
            An instance of LinkList created from the CIF DataContainer.
        """
        L = []
        cif_category = dc.getObj(Link._CIF_CategoryName)
        if cif_category is None:
            return cls(L)  # Return empty LinkList if category not found
        for i in range(len(cif_category)):
            for key, valset in Link._CIF_CategoryElementTypes.items():
                objTypeid = cif_category.getValue(key, i)
                if objTypeid in valset:
                    this_link = Link(CIFdict(cif_category, i))
                    L.append(this_link)
        return cls(L)

    def describe(self) -> str:
        """
        Returns a string description of the LinkList.
        
        Returns
        -------
        str
            A string describing the number of links in the list.
        """
        return f'<LinkList with {len(self)} links>'

    def assign_residues(self, Residues: "ResidueList") -> tuple["ResidueList", "LinkList"]:
        """
        Assigns residue and atom pointers to each link; sets up the up and down links of both
        residues so that linked residue objects can reference one another; flags residues from
        list of residues passed in that are not assigned to any links

        Parameters
        ----------
        Residues: ResidueList
            list of residues to assign to links; this is typically a list of all residues in the
            structure, but it can also be a list of residues that are not linked to any other
            residues, such as when reading in a set of links from a pre-built psf file.

        Returns
        -------
        Residues: ResidueList
            list of residues from Residues that are not used for any assignments
        """
        logger.debug(f'Links: Assigning residues from list of {len(Residues)} residues')
        ignored_by_ptnr1 = self.assign_objs_to_attr('residue1', Residues, chainID='chainID1', resid='resid1')
        ignored_by_ptnr2 = self.assign_objs_to_attr('residue2', Residues, chainID='chainID2', resid='resid2')
        for link in self.data:
            link.segtype1 = link.residue1.segtype
            link.segtype2 = link.residue2.segtype
            # shortcodes don't provide resnames, so set them here
            if link.residue1 is not None and link.resname1 is None:
                link.resname1 = link.residue1.resname
            if link.residue2 is not None and link.resname2 is None:
                link.resname2 = link.residue2.resname
            try:
                link.residue1.link_to(link.residue2, link)
            except:
                raise ValueError(f'Bad residue in link')
            if link.patchname is not None and len(link.patchname) > 0:
            # this link was most likely created when reading in patch records from a set of REMARKS in a pre-built psf file.
            # we need to get the precise atom names for this patch
                continue

            link.atom1 = link.residue1.atoms.get(lambda x: x.name == link.name1 and x.altloc == link.altloc1)
            link.atom2 = link.residue2.atoms.get(lambda x: x.name == link.name2 and x.altloc == link.altloc2)
            link.segtype1 = link.residue1.segtype
            link.segtype2 = link.residue2.segtype
            # shortcodes don't provide resnames, so set them here
            if link.residue1 is not None and link.resname1 is None:
                link.resname1 = link.residue1.resname
            if link.residue2 is not None and link.resname2 is None:
                link.resname2 = link.residue2.resname
            link.set_patchname()
        # do cross-assignment to find true orphan links and dangling links
        orphan_1 = ignored_by_ptnr1.assign_objs_to_attr('residue2', Residues, chainID='chainID2', resid='resid2')
        orphan_2 = ignored_by_ptnr2.assign_objs_to_attr('residue1', Residues, chainID='chainID1', resid='resid1')
        orphans = orphan_1 + orphan_2
        rlist = []
        for link in ignored_by_ptnr1:
            rlist, list = link.residue2.get_down_group()
            rlist.insert(0, link.residue2)
            for r in rlist:
                Residues.remove(r)
        return Residues.__class__(rlist), self.__class__(ignored_by_ptnr1 + ignored_by_ptnr2)

    def remove_links_to(self, r) -> 'LinkList':
        """
        Remove all links to a specific residue.

        Parameters
        ----------
        Residue: Residue
            The residue to remove links to.

        Returns
        -------
        LinkList
            A new LinkList containing the removed links.
        """
        removed_links = self.__class__([])
        for link in self.data:
            if link.residue1 == r or link.residue2 == r:
                removed_links.append(link)
        for badlink in removed_links:
            self.remove(badlink)
        logger.debug(f'Removed {len(removed_links)} links to residue {str(r)}')
        return removed_links

    # def prune_mutations(self, Mutations: 'MutationList', Segments: 'SegmentList'):
    #     """
    #     Prune off any links and associated objects as a result of mutations

    #     Parameters
    #     ----------
    #     Mutations: MutationList
    #         list of mutations to prune off links; these are typically the mutations that were applied
    #         to the structure before the links were created, so they are not part of the original
    #         structure, but they are part of the current structure.
    #     Segments: SegmentList
    #         Current list of segments in the structure; it is from this list that segments are removed if they become empty as a result of pruning off links that are pruned off due to mutations.

    #     Returns
    #     -------
    #     pruned: dict
    #         A dictionary containing lists of pruned residues, links, and segments.
        
    #         - ``residues``: list of Residue objects that were pruned
    #         - ``links``: list of Link objects that were pruned
    #         - ``segments``: list of Segment objects that were pruned
    #     """
    #     pruned = {'residues': [], 'links': self.__class__([]), 'segments': []}
    #     for m in Mutations.data:
    #         left = self.get(lambda x: x.chainID1 == m.chainID and x.resid1 == m.resid) # links for which partner 1 is the mutation
    #         right = self.get(lambda x: x.chainID2 == m.chainID and x.resid2 == m.resid) # links for which partner 2 is the mutation
    #         if left:  # this is a link in which this mutation is the partner 1
    #             self.remove(left) # get rid of this link
    #             # we need to remove residue2 and everything downstream
    #             # remove downstream residues!
    #             rlist, llist = left.residue2.get_down_group()
    #             rlist.insert(0, left.residue2)
    #         elif right: # this is a link in which this mutation is the right member (should be very rare)
    #             self.remove(right)
    #             rlist, llist = right.residue2.get_down_group()
    #             rlist.insert(0, right.residue2)
    #         if rlist and llist:
    #             # logger.debug(f'Deleting residues down from and including {str(rlist[0])} due to a mutation')
    #             S = Segments.get_segment_of_residue(rlist[0])
    #             for r in rlist:
    #                 # logger.debug(f'...{str(r)}')
    #                 S.residues.remove(r)
    #                 pruned['residues'].append(r)
    #             if len(S.residues) == 0:
    #                 # logger.debug(f'All residues of {S.psfgen_segname} are deleted; {S.psfgen_segname} is deleted')
    #                 Segments.remove(S)
    #                 pruned['segments'].append(S)
    #             for l in llist:
    #                 self.data.remove(l)
    #                 pruned['links'].append(l)
    #     return pruned

    def apply_segtypes(self, map):
        """
        Apply segtype values to each of the two residues using the map
        
        Parameters
        ----------
        map: dict
            map of segtypes for given resnames
        """
        self.map_attr('segtype1','resname1', map)
        self.map_attr('segtype2','resname2', map)

    def report(self) -> str:
        """
        Report the string representation of each link in the list.
        """
        return "\n".join([str(l) for l in self])
    
def ic_reference_closest(res12: list["Residue"], ICmaps: list[dict]) -> str:
    """
    Given the two Residues in res12 and the maps in ICmaps, 
    return the mapping key to which the given IC values are
    closest in a Euclidean sense.

    This method will identify the four atoms of the IC and reference
    them directly when calling the measure_dihedral function. 
    
    The list of computed dihedral values from the set of atoms is a "point"
    in "IC-space", and each patch has its own "reference point" in this space.

    The reference point to which the point is closest is identified as the 
    desired result.

    Parameters
    ----------
    res12: list
       exactly two Residue objects which must have lists of atoms attributes
    ICMaps : list of dict
        A list of dictionaries, each with the following structure:

        +--------------+-----------------------------------------------------------+
        | Key          | Description                                               |
        +==============+===========================================================+
        | ICatomnames  | List of 4 atom names as they appear in the CHARMM FF IC.  |
        +--------------+-----------------------------------------------------------+
        | mapping      | Dictionary mapping patch names to IC values.              |
        +--------------+-----------------------------------------------------------+
    """
    for ic in ICmaps:
        # logger.debug(f'icmap {ic}')
        ic['atoms'] = []
        for n in ic['ICatomnames']:
            r = int(n[0])-1
            an = n[1:]
            at = res12[r].atoms.get(lambda x: x.name == an)
            ic['atoms'].append(at)
            # logger.debug(f'Assigned atom {at.name} of {at.resname}{at.resseqnum}')
    map_points = {}
    the_point = []
    for ic in ICmaps:
        value = measure_dihedral(*(ic['atoms'])) * 180.0 / np.pi
        # logger.debug(f'{ic["ICatomnames"]} value {value:.2f}')
        the_point.append(value)
        for m, v in ic['mapping'].items():
            if not m in map_points:
                map_points[m] = []
            map_points[m].append(v)
    for k, v in map_points.items():
        map_points[k] = np.array(v)
    the_point = np.array(the_point)
    # logger.debug(f'ic the point: {the_point}')
    # calculate Euclidean distance adhering to the periodicity
    # of dihedral-angle space
    displacements = {k: (the_point - v) for k, v in map_points.items()}
    for n, d in displacements.items():
        for i in range(len(d)):
            if d[i] < 180.0:
                d[i] += 180.0
            if d[i] > 180.0:
                d[i] -= 180.0
    norms = {k: np.linalg.norm(d) for k, d in displacements.items()}
    # logger.debug(f'norms {norms}')
    the_one = [k for (k, v) in sorted(norms.items(), key=lambda x: x[1])][0]
    # logger.debug(f'returning {the_one}')
    return the_one

