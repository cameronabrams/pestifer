# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
Defines the ResiduePlaceholder and Residue classes for handling residues in a molecular structure.
"""

from __future__ import annotations

import logging

import numpy as np

from argparse import Namespace

from pydantic import Field
from typing import ClassVar, Callable

from pidibble.baserecord import BaseRecord
from pidibble.pdbrecord import PDBRecord, PDBRecordDict
from mmcif.api.PdbxContainers import DataContainer

from pestifer.util.stringthings import parse_filter_expression

from .atom import Atom, AtomList, Hetatm
from ..objs.link import Link
from ..objs.link import LinkList

from ..core.baseobj import BaseObj, BaseObjList
from ..core.labels import Labels
from ..objs.deletion import DeletionList
from ..objs.insertion import InsertionList
from ..objs.patch import Patch
from ..objs.seqadv import Seqadv, SeqadvList
from ..objs.ssbond import SSBond
from ..objs.substitution import SubstitutionList
from ..util.cifutil import CIFdict
from ..util.coord import positionN
from .stateinterval import StateInterval, StateIntervalList

from ..objs.resid import ResID

logger = logging.getLogger(__name__)

class ResiduePlaceholder(BaseObj):
    """
    A class for handling missing residues in a molecular structure.
    This class represents residues that are not present in the structure, such as those that are missing
    due to low resolution or other reasons. It is used to track residues that are expected to be present
    but are not resolved in the coordinate file.
    """

    _required_fields = {'resname', 'resid', 'chainID', 'resolved', 'segname', 'segtype'}
    """
    Required attributes for ResiduePlaceholder.
    
    Attributes
    ----------
    resname : str
        The residue name.
    resid : ResID
        The residue ID (sequence number with optional 1-character insertion code).
    chainID : str
        The chain ID.
    resolved : bool
        Indicates whether the residue is resolved (True) or not (False).
    segname : str
        Segment name of the residue.
    segtype : str
        The segment type.
    """

    _optional_fields = {'model', 'obj_id', 'auth_asym_id', 'auth_comp_id', 'pdb_ins_code',
                        'auth_seq_id', 'empty', 'recordname', 'asym_chainID', 'ORIGINAL_ATTRIBUTES'}
    """
    Optional attributes for ResiduePlaceholder.
    
    Attributes
    ----------
    model : int
        The model number, if applicable.
    obj_id : str
        The residue ID.
    auth_asym_id : str
        The author asymmetry ID, if applicable.
    auth_comp_id : str
        The author component ID, if applicable.
    auth_seq_id : int
        The author sequence ID, if applicable.
    empty : bool
        Indicates whether the residue is empty (True) or not (False).
    recordname : str
        The PDB record name for the residue.
    asym_chainID : str
        The chain ID of this residue's image in the asymmetric unit, if applicable.
    """

    _ignore_fields = {'recordname', 'empty'}
    """
    Fields that are ignored when comparing ResiduePlaceholder objects.
    """
    
    resname: str = Field(..., description="The residue name.")
    resid: ResID = Field(..., description="The residue ID.")
    chainID: str = Field(..., description="The chain ID.")
    resolved: bool = Field(..., description="Indicates whether the residue is resolved (True) or not (False).")
    segtype: str = Field(..., description="The segment type.")
    segname: str = Field(..., description="The segment name.")
    model: int = Field(1, description="The model number, if applicable.")
    obj_id: int | None = Field(None, description="The object id, if applicable.")
    auth_asym_id: str | None = Field(None, description="The author asymmetry ID, if applicable.")
    auth_comp_id: str | None = Field(None, description="The author component ID, if applicable.")
    auth_seq_id: int | None = Field(None, description="The author sequence ID, if applicable.")
    pdb_ins_code: str | None = Field(None, description="The PDB insertion code, if applicable.")

    empty: bool = Field(False, description="Indicates whether the residue is empty (True) or not (False).")
    recordname: str = Field('REMARK.465', description="The PDB record name for the residue.")

    asym_chainID: str | None = Field(None, description="The original chain ID in the asymmetric unit, if applicable.")
    ORIGINAL_ATTRIBUTES: dict = Field(default_factory=dict, description="Stash of original attributes for the residue, if it is renamed.")

    _yaml_header: ClassVar[str] = 'missings'
    """
    YAML header for ResiduePlaceholder.
    """

    _PDB_keyword: ClassVar[str] = 'REMARK.465'
    """
    PDB keyword for ResiduePlaceholder.
    """
    _PDB_table_of_keyword: ClassVar[str] = 'MISSING'
    """
    PDB table keyword for ResiduePlaceholder.
    """

    _CIF_CategoryName: ClassVar[str] = 'pdbx_unobs_or_zero_occ_residues'
    """
    mmCIF category name for ResiduePlaceholder.
    """

    def shortcode(self):
        return f"{self.model}:{self.chainID}:{self.resname}{self.resid.resid}"

    @staticmethod
    def _from_shortcode(shortcode: str) -> dict:
            # i:C:RD; e.g., # 1:A:ALA123 or 1:A:ALA123A or A:A123
            # i model number, optional, default 1
            # C chain ID
            # R one or three-letter residue code
            # D is resid
            tokens = shortcode.split(':')
            has_modelnum = len(tokens) > 2
            model_idx = -1
            chain_idx = -1
            model = 1
            chainID = 'A'
            if has_modelnum:
                model_idx = 0
                model = int(tokens[model_idx])
            has_chainID = len(tokens) > 1
            if has_chainID:
                chain_idx = model_idx + 1
                chainID = tokens[chain_idx]
            res_idx = chain_idx + 1
            rescode_and_ID = tokens[res_idx]
            # rescode are the first alphabetic characters of rescode_and_ID.
            rescode_attempt = rescode_and_ID[:3].upper()
            # if this is a single-letter code, expect the second character to be a digit
            if rescode_attempt[1].isdigit():
                rescode = rescode_attempt[0]
                resid = ResID(rescode_and_ID[1:])
                resname = Labels.res_123.get(rescode, rescode_attempt)
            else:
                resname = rescode_attempt
                resid = ResID(rescode_and_ID[3:])
            input_dict = {
                'model': model,
                'resname': resname,
                'chainID': chainID,
                'asym_chainID': chainID,
                'segname' : 'UNSET',
                'resid': resid,
                'resolved': False,
                'segtype': 'UNSET',
            }
            return input_dict

    @staticmethod
    def _from_pdbrecord(pdbrecord: PDBRecord | BaseRecord) -> dict:
        """
        Converts a PDBRecord to a dictionary of parameters for ResiduePlaceholder.
        
        Parameters
        ----------
        pdbrecord : PDBRecord
            A PDBRecord instance containing residue information.
        
        Returns
        -------
        dict
            A dictionary containing the parameters for ResiduePlaceholder.
        """
        if pdbrecord.modelNum is None:
            pdbrecord.modelNum = 1
        elif type(pdbrecord.modelNum) == str and pdbrecord.modelNum.isdigit():
            pdbrecord.modelNum = int(pdbrecord.modelNum)
        elif type(pdbrecord.modelNum) == str and len(pdbrecord.modelNum) == 0:
            pdbrecord.modelNum = 1
        input_dict = {
            # pidbble/resources/pdb_format.yaml line 846 to 850
            'model': pdbrecord.modelNum,
            'resname': pdbrecord.resName,
            'chainID': pdbrecord.chainID,
            'asym_chainID' : pdbrecord.chainID,
            'resid': ResID(pdbrecord.seqNum, pdbrecord.iCode),
            'resolved': False,
            'segtype': 'UNSET',
            'segname': 'UNSET'
        }
        return input_dict

    @staticmethod
    def _from_cifdict(cifdict: CIFdict) -> dict:
        mn = cifdict['pdb_model_num']
        if type(mn) == str and mn.isdigit():
            nmn = int(mn)
        else:
            nmn = 1
        input_dict = {
            'model': nmn,
            'resname': cifdict['label_comp_id'],
            'chainID': cifdict['label_asym_id'],
            'asym_chainID': cifdict['label_asym_id'],
            'resid': ResID(int(cifdict['label_seq_id'])),
            'resolved': False,
            'segtype': 'UNSET',
            'segname': 'UNSET',
            'auth_asym_id': cifdict['auth_asym_id'],
            'auth_comp_id': cifdict['auth_comp_id'],
            'auth_seq_id': int(cifdict['auth_seq_id']),
            'pdb_ins_code': cifdict['pdb_ins_code'],
        }
        return input_dict

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        if args and isinstance(args[0], str):
            return cls._from_shortcode(args[0])
        elif args and isinstance(args[0], PDBRecord | BaseRecord):
            return cls._from_pdbrecord(args[0])
        elif args and isinstance(args[0], CIFdict):
            return cls._from_cifdict(args[0])
        return super()._adapt(*args, **kwargs)

    def __str__(self):
        return f'{self.chainID}_{self.resname}{self.resid.resid}*'

    def pdb_line(self):
        """
        Returns a PDB line representation of the :class:`ResiduePlaceholder`.
        The line is formatted according to the PDB standard for missing residues.
        
        Returns
        -------
        str
            A string representing the :class:`ResiduePlaceholder` in PDB format.
        """
        record_name, code = ResiduePlaceholder._PDB_keyword.split('.')
        return '{:6s}{:>4d}   {:1s} {:3s} {:1s} {:>6s}'.format(record_name,
        code, self.model, self.resname, self.chainID, self.resid.pdbresid)

class ResiduePlaceholderList(BaseObjList[ResiduePlaceholder]):
    """
    A class for handling lists of :class:`ResiduePlaceholder` objects.
    This class is used to manage collections of residues that are not present in the molecular structure,
    such as missing residues in a PDB file.

    This class does not add anything beyond the :class:`~pestifer.core.baseobj.BaseObjList` class.
    """
    def describe(self):
        return f'<ResiduePlaceholderList with {len(self)} residues>'
    
    @classmethod
    def from_pdb(cls, parsed: PDBRecordDict) -> "ResiduePlaceholderList":
        """
        Create an ResiduePlaceholderList from parsed PDB data.
        
        Parameters
        ----------
        parsed : PDBRecordDict
            A dictionary containing parsed PDB records.
        
        Returns
        -------
        ResiduePlaceholderList
            A new instance of ResiduePlaceholderList containing the missing residues.
        """
        if ResiduePlaceholder._PDB_keyword not in parsed:
            return cls([])
        return cls([ResiduePlaceholder(p) for p in parsed[ResiduePlaceholder._PDB_keyword].tables[ResiduePlaceholder._PDB_table_of_keyword]])
    
    @classmethod
    def from_cif(cls, cif_data: DataContainer) -> "ResiduePlaceholderList":
        """
        Create an ResiduePlaceholderList from CIF data.
        
        Parameters
        ----------
        cif_data : DataContainer
            A DataContainer instance containing CIF data.
        
        Returns
        -------
        ResiduePlaceholderList
            A new instance of ResiduePlaceholderList containing the missing residues.
        """
        obj = cif_data.getObj(ResiduePlaceholder._CIF_CategoryName)
        if obj is None:
            return cls([])
        return cls([ResiduePlaceholder(CIFdict(obj, i)) for i in range(len(obj))])

    def apply_inclusion_logics(self, inclusion_logics: list[str] = []):
        ignored_missing_residue_count = 0
        if len(inclusion_logics) == 0:
            return 0
        keep_missing_residues = ResiduePlaceholderList([])
        for expression in inclusion_logics:
            logger.debug(f'Applying inclusion logic: {expression}')
            filter_func = parse_filter_expression(expression)
            keep_missing_residues.extend(list(filter(filter_func, self.data)))
        logger.debug(f'Keeping {len(keep_missing_residues)} missing residues.')
        num_kept_missing_residues = len(keep_missing_residues)
        if num_kept_missing_residues > 0:
            self.data = keep_missing_residues
        return ignored_missing_residue_count

    def apply_exclusion_logics(self, exclusion_logics: list[str] = []) -> int:
        if len(exclusion_logics) == 0:
            return 0
        all_ignored_missing_residues = ResiduePlaceholderList([])
        for expression in exclusion_logics:
            logger.debug(f'Excluding missing residues with expression: {expression}')
            filter_func = parse_filter_expression(expression)
            ignored_missing_residues = list(filter(filter_func, self.data))
            all_ignored_missing_residues.extend(ignored_missing_residues)
        logger.debug(f'Removing {len(all_ignored_missing_residues)} ignored missing residues.')
        for residue in all_ignored_missing_residues:
            self.remove_instance(residue)
        return len(all_ignored_missing_residues)
    
class Residue(ResiduePlaceholder):
    """
    A class for handling residues in a molecular structure.
    This class extends the :class:`ResiduePlaceholder` class to include additional functionality for managing residues.
    """

    _required_fields = ResiduePlaceholder._required_fields | {'atoms'}
    """
    Required attributes for :class:`~pestifer.molecule.residue.Residue` in addition to those defined for :class:`~pestifer.molecule.residue.ResiduePlaceholder`.

    Attributes
    ----------
    atoms : :class:`~pestifer.molecule.atom.AtomList`
        A list of :class:`~pestifer.molecule.atom.Atom` objects that belong to this residue.
    """
    atoms: AtomList = Field(default_factory=AtomList, description="A list of atoms that belong to this residue.")

    # _optional_fields = ResiduePlaceholder._optional_fields | {'uplink', 'downlink'}
    _optional_fields = ResiduePlaceholder._optional_fields | {'up', 'down', 'uplink', 'downlink'}
    """
    Optional attributes for :class:`~pestifer.molecule.residue.Residue` in addition to those defined for :class:`~pestifer.molecule.residue.ResiduePlaceholder`.

    - ``up``: list
        A list of residues that are linked to this residue in an upstream direction.
    - ``down``: list
        A list of residues that are linked to this residue in a downstream direction.
    - ``uplink``: list
        A list of links to residues that are connected to this residue in an upstream direction.
    - ``downlink``: list
        A list of links to residues that are connected to this residue in a downstream direction.
    """
    up: 'ResidueList' = Field(default_factory=lambda: ResidueList(), description="A list of residues linked to this residue in an upstream direction.")
    down: 'ResidueList' = Field(default_factory=lambda: ResidueList(), description="A list of residues linked to this residue in a downstream direction.")
    uplink: LinkList = Field(default_factory=LinkList, description="A list of links to residues connected to this residue in an upstream direction.")
    downlink: LinkList = Field(default_factory=LinkList, description="A list of links to residues connected to this residue in a downstream direction.")

    _ignore_fields = ResiduePlaceholder._ignore_fields | {'atoms', 'up', 'down', 'uplink', 'downlink'}
    """
    Attributes to ignore when comparing Residue objects.
    This includes the attributes defined in :class:`~pestifer.molecule.residue.ResiduePlaceholder` as well as the additional attributes defined for :class:`~pestifer.molecule.residue.Residue`.

    - ``atoms``: :class:`~pestifer.molecule.atom.AtomList`
        A list of :class:`~pestifer.molecule.atom.Atom` objects that belong to this residue.
    - ``up``: list
        A list of residues that are linked to this residue in an upstream direction.
    - ``down``: list
        A list of residues that are linked to this residue in a downstream direction.
    - ``uplink``: list
        A list of links to residues that are connected to this residue in an upstream direction.
    - ``downlink``: list
        A list of links to residues that are connected to this residue in a downstream direction.
    """

    _PDB_keyword: ClassVar[str] = None
    _CIF_CategoryName: ClassVar[str] = None
    _yaml_header: ClassVar[str] = None

    _counter: ClassVar[int] = 0

    _chainIDmap_label_to_auth = {}

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Converts various input types into a dictionary of parameters for Residue.
        
        Parameters
        ----------
        args : tuple
            A tuple containing the input data to be converted.
        
        Returns
        -------
        dict
            A dictionary containing the parameters for Residue.
        
        Raises
        ------
        TypeError
            If the input type is not supported.
        """
        if args and isinstance(args[0], str):
            input_dict = ResiduePlaceholder._adapt(args[0])
            input_dict['recordname'] = 'UNSET'
            input_dict['atoms'] = AtomList()
            return input_dict
        elif args and isinstance(args[0], ResiduePlaceholder):
            input_dict = {'resname': args[0].resname,
                          'resid': ResID(args[0].resid),
                          'chainID': args[0].chainID,
                          'asym_chainID': args[0].asym_chainID,
                          'segtype': args[0].segtype,
                          'segname': args[0].segname,
                          'model': args[0].model,
                          'obj_id': args[0].obj_id,
                          'auth_asym_id': args[0].auth_asym_id,
                          'auth_comp_id': args[0].auth_comp_id,
                          'auth_seq_id': args[0].auth_seq_id,
                          'atoms': AtomList(),
                          'resolved': False,
                          'recordname': 'UNSET'}
            return input_dict
        elif args and isinstance(args[0], PDBRecord | BaseRecord):
            input_dict = ResiduePlaceholder._adapt(args[0])
            input_dict['atoms'] = AtomList()
            input_dict['resolved'] = False
            input_dict['recordname'] = 'UNSET'
            input_dict['segname'] = input_dict.get('chainID','UNSET')
            return input_dict
        elif args and isinstance(args[0], CIFdict):
            input_dict = ResiduePlaceholder._adapt(args[0])
            input_dict['atoms'] = AtomList()
            input_dict['resolved'] = False
            input_dict['recordname'] = 'UNSET'
            input_dict['segname'] = input_dict.get('chainID','UNSET')
            return input_dict
        elif args and isinstance(args[0], Atom | Hetatm):
            a = args[0]
            # logger.debug(f'Adapting Atom {str(a)} to Residue (resid={a.resid})')
            input_dict = dict(resname=a.resname,
                              resid=a.resid, chainID=a.chainID,
                              asym_chainID=a.chainID,
                              atoms=AtomList([a]),
                              segtype=Labels.segtype_of_resname[a.resname],
                              segname=a.segname,
                              auth_asym_id=getattr(a, 'auth_asym_id', None),
                              auth_comp_id=getattr(a, 'auth_comp_id', None),
                              auth_seq_id=getattr(a, 'auth_seq_id', None),
                              resolved=True)
            input_dict['recordname'] = 'UNSET'
            return input_dict
        return super()._adapt(*args, **kwargs)
            
    def __str__(self):
        return f'{self.chainID}_{self.resname}{self.resid.resid}'

    def __lt__(self, other):
        """
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "less than" the other; i.e., N-terminal to.  This is used for sorting residues in a sequence.
        The comparison is based on the residue sequence number and insertion code.
        If the other object is a string, it is split into residue sequence number and insertion code
        using the :func:`~pestifer.core.stringthings.split_ri` function.

        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`, str
            The other residue or string to compare with.
        
        Returns
        -------
        bool
            True if this residue is less than the other residue, False otherwise.
        """
        if isinstance(other, Residue):
            return self.resid < other.resid
        elif isinstance(other, str):
            return self.resid < ResID(other)
        else:
            raise TypeError(f"Cannot compare Residue with {type(other)}")
    
    def __gt__(self, other):
        """
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "greater than" the other; i.e., C-terminal to.  This is used for sorting residues in a sequence. If the other object is a string, it is split into residue sequence number and insertion code using the :func:`~pestifer.core.stringthings.split_ri` function.

        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`, str
            The other residue or string to compare with.

        Returns
        -------
        bool
            True if this residue is greater than the other residue, False otherwise.
        """
        if isinstance(other, Residue):
            return self.resid > other.resid
        elif isinstance(other, str):
            return self.resid > ResID(other)
        else:
            raise TypeError(f"Cannot compare Residue with {type(other)}")
        
    def __eq__(self, other):
        """
        Check if this residue is equal to another residue or a string representation of a residue.
        This is used to determine if two residues are the same in terms of their sequence position.

        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`, str
            The other residue or string to compare with.

        Returns
        -------
        bool
            True if this residue has the same residue sequence number and insertion code as the other residue,
            False otherwise.
        """
        if isinstance(other, Residue):
            # comparing to another residue requires same resid and chainID
            return self.resid == other.resid and self.chainID == other.chainID
        elif isinstance(other, ResID):
            # comparing just to a Resid requires just same resid
            return self.resid == other
        elif isinstance(other, str):
            # str interpreted as resid
            return self.resid == ResID(other)
        else:
            raise TypeError(f"Cannot compare Residue with {type(other)}")

    def __ne__(self, other):
        """
        Check if this residue is not equal to another residue or a string representation of a residue.
        This is used to determine if two residues are different in terms of their sequence position.

        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`, str
            The other residue or string to compare with.

        Returns
        -------
        bool
            True if this residue does not have the same residue sequence number and insertion code as the other residue,
            False otherwise.
        """
        return self.resid != other.resid if isinstance(other, Residue) else self.resid != ResID(other)

    def __le__(self,other):
        """
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "less than or equal to" the other; i.e., N-terminal to or same as.  This is used for sorting residues in a sequence.

        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`, str
            The other residue or string to compare with.

        Returns
        -------
        bool
            True if this residue is less than or equal to the other residue, False otherwise.
        """
        if isinstance(other, Residue):
            return self.resid == other.resid or self.resid < other.resid
        elif isinstance(other, str):
            return self.resid == ResID(other) or self.resid < ResID(other)
        raise TypeError(f"Cannot compare Residue with {type(other)}")
    
    def __ge__(self,other):
        """
        Compare this residue with another residue or a string representation of a residue to determine if this residue is "greater than or equal to" the other; i.e., C-terminal to or same as.  This is used for sorting residues in a sequence.

        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`, str
            The other residue or string to compare with.

        Returns
        -------
        bool
            True if this residue is greater than or equal to the other residue, False otherwise.
        """
        if isinstance(other, Residue):
            return self.resid == other.resid or self.resid > other.resid
        elif isinstance(other, str):
            return self.resid == ResID(other) or self.resid > ResID(other)
        raise TypeError(f"Cannot compare Residue with {type(other)}")


    def same_resid(self,other):
        """
        Check if this residue has the same residue sequence number and insertion code as another residue or a
        string representation of a residue. This is used to determine if two residues are the same in terms of their sequence position.
        
        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`, str
            The other residue or string to compare with.
            
        Returns
        -------
        bool
            True if this residue has the same residue sequence number and insertion code as the other residue,
            False otherwise.
        """
        if type(other) == type(self):
            return self.resid == other.resid
        elif type(other) == str:
            return self.resid == ResID(other)
        return False

    def add_atom(self, a: Atom):
        """
        Add an atom to this residue if it matches the residue's resid, residue name,
        and chain ID. This method is used to build a residue from its constituent
        atoms, ensuring that all atoms in the residue share the same resid, residue name,
        and chain ID.

        Parameters
        ----------
        a : :class:`~pestifer.molecule.atom.Atom`
            The atom to be added to the residue.
        """
        # logger.debug(f'resid: {self.resid} == {a.resid}; resname: {self.resname} == {a.resname}; chainID: {self.chainID} == {a.chainID}')
        if self.resid == a.resid and self.resname == a.resname and self.chainID == a.chainID:
            self.atoms.append(a)
            # logger.debug(f'Atom {repr(a)} added to residue {repr(self)}')
            return True
        # logger.debug(f'Atom {repr(a)} not added to residue {repr(self)}')
        return False

    def set_chainID(self, chainID):
        """
        Set the chain ID for this residue and all its constituent atoms.
        This method updates the chain ID of the residue and all atoms within it to ensure consistency
        across the residue's structure. 
        
        Parameters
        ----------
        chainID : str
            The new chain ID to be set for the residue and its atoms.
        """
        self.set(asym_chainID=self.chainID)
        self.set(chainID=chainID)
        for a in self.atoms.data:
            a.set(chainID=chainID)
    
    def set_segname(self, segname: str):
        self.set(segname=segname)

    def link_to(self, other: Residue, link: Link):
        """
        Link this residue to another residue in a downstream direction.
        This method establishes a connection between this residue and another residue, allowing for
        traversal of the molecular structure in a downstream direction. It also updates the upstream
        and downstream links for both residues to maintain the connectivity information.
        
        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`
            The other residue to link to.
        link : :class:`~pestifer.objs.link.Link`
            A link object representing the connection between the two residues.
        """
        assert type(other)==type(self),f'type of other is {type(other)}; expected {type(self)}'
        self.down.append(other)
        self.downlink.append(link)
        other.up.append(self)
        other.uplink.append(link)

    def unlink(self, other: Residue, link: Link):
        """
        Unlink this residue from another residue in a downstream direction.
        This method removes the connection between this residue and another residue, effectively
        severing the link between them. It updates the upstream and downstream links for both residues
        to reflect the removal of the connection.
        
        Parameters
        ----------
        other : :class:`~pestifer.molecule.residue.Residue`
            The other residue to unlink from.
        link : :class:`~pestifer.objs.link.Link`
            A link object representing the connection between the two residues.
        """
        self.down.remove(other)
        self.downlink.remove(link)
        other.up.remove(self)
        other.uplink.remove(link)

    def ri(self):
        """
        Returns the residue sequence number and insertion code as a string.
        """
        return self.resid.resid
    
    def get_down_group(self):
        """
        Get the downstream group of residues.
        This method traverses the downstream links of this residue and collects all residues
        in the downstream direction.
        """
        res=ResidueList()
        lin=LinkList()
        for d,dl in zip(self.down,self.downlink):
            res.append(d)
            lin.append(dl)
            # logger.debug(f'{str(self)}->{str(d)}')
            tres,tlin=d.get_down_group()
            res.extend(tres)
            lin.extend(tlin)
        return res,lin

class ResidueList(BaseObjList[Residue]):
    """
    A class for handling lists of :class:`~pestifer.molecule.residue.Residue` objects.
    This class extends the :class:`~pestifer.core.baseobj.AncestorAwareObjList` to manage collections of residues in a molecular structure.
    It provides methods for initializing the list from various input types, indexing residues, and performing operations
    such as mapping chain IDs, retrieving residues and atoms, and handling residue ranges.
    """
    def describe(self):
        """
        Returns a string description of the ResidueList object, including the number of residues it contains.
        
        Returns
        -------
        str
            A string representation of the ResidueList object.
        """
        return f'<ResidueList with {len(self)} residues>'
    
    @classmethod
    def from_residuegrouped_atomlist(cls, atoms: AtomList):
        """
        Construct a list of residues from a list of Atoms ordered and grouped into residues
        """
        R = []
        for atom in atoms.data:
            current_residue = None if len(R) == 0 else R[-1]
            if not current_residue or not current_residue.add_atom(atom):
                R.append(Residue(atom))
                # logger.debug(f'Created new residue {R[-1].resname}_{R[-1].resid.resid} with first atom {atom.name}')
            # logger.debug(f'len(R) {len(R)} len(atoms) {len(atoms)}')

        # for r in R:
        #     logger.debug(f'Adding residue {r.resname}{r.resid.resid} with {len(r.atoms)} atoms')
        return cls(R)
    
    @classmethod
    def from_ResiduePlaceholderlist(cls, input_list: ResiduePlaceholderList):
        """
        Create a ResidueList from an ResiduePlaceholderList.
        
        Parameters
        ----------
        input_list : ResiduePlaceholderList
            A list of ResiduePlaceholder objects to convert into a ResidueList.
        
        Returns
        -------
        ResidueList
            An instance of ResidueList initialized with the residues created from the empty residues.
        """
        R = [Residue(m) for m in input_list.data]
        return cls(R)

    def map_chainIDs_label_to_auth(self):
        """
        Create a mapping from chain IDs in the label (e.g., PDB format) to the author chain IDs.
        """
        self._chainIDmap_label_to_auth = {}
        for r in self.data:
            if hasattr(r,'auth_asym_id'):
                label_Cid = r.chainID
                auth_Cid = r.auth_asym_id
                if not label_Cid in self._chainIDmap_label_to_auth:
                    self._chainIDmap_label_to_auth[label_Cid] = auth_Cid
    
    def atom_serials(self, as_type=str):
        """
        Get a list of atom serial numbers for all atoms in the residues.
        
        Parameters
        ----------
        as_type : type, optional
            The type to which the serial numbers should be converted. Default is str.

        Returns
        -------
        list
            A list of atom serial numbers.
        """
        serlist=[]
        for res in self.data:
            for a in res.atoms:
                serlist.append(as_type(a.serial))
        return serlist

    def atom_resseqnums(self, as_type=str):
        """
        Get a list of residue sequence numbers for all atoms in the residues.
        
        Parameters
        ----------
        as_type : type, optional
            The type to which the residue sequence numbers should be converted. Default is str.
        
        Returns
        -------
        list
            A list of residue sequence numbers.
        """
        rlist=[]
        for res in self.data:
            for a in res.atoms:
                rlist.append(as_type(a.resseqnum))
        return rlist
    
    def atom_resids(self, as_type=ResID):
        """
        Get a list of residue IDs for all atoms in the residues.
        
        Parameters
        ----------
        as_type : type, optional
            The type to which the residue IDs should be converted. Default is str.
        
        Returns
        -------
        list
            A list of residue IDs.
        """
        rlist=[]
        for res in self.data:
            for a in res.atoms:
                rlist.append(as_type(a.resid))
        return rlist

    def caco_str(self, upstream_reslist: ResidueList, seglabel: str, tmat: np.ndarray) -> str:
        """
        Generate a string representation of the ``caco`` command to position
        the N atom of a model-built residue based on the coordinates of the previous residue.  Uses the :func:`~pestifer.util.coord.positionN` function to calculate the position of the N atom based on the upstream residue's coordinates and the transformation matrix.
        
        Parameters
        ----------
        upstream_reslist : list
            A list of residues that are upstream of the current residue.
        seglabel : str
            The segment label for the residue.
        molid_varname : str
            The variable name for the molecule ID. (unused for now)
        tmat : :class:`numpy.ndarray`
            A transformation matrix to apply to the coordinates of the residue.

        Returns
        -------
        str
            A string representing the VMD/psfgen ``caco`` command for positioning the residue.
        """
        r0 = self[0]
        ur = upstream_reslist[-1]
        rN = positionN(ur, tmat)
        logger.debug(f'caco {rN}')
        return f'coord {seglabel} {r0.resid.resid} N {{{rN[0]:.5f} {rN[1]:.5f} {rN[2]:.5f}}}'
    
    def resrange(self, rngrec):
        """
        Yield residues from specified residue range.
        
        Parameters
        ----------
        rngrec : :class:`argparse.Namespace`
            A record defining the range of residues to retrieve. It must have the following attributes:

            - `chainID`: The chain ID of the residues.
            - `resid1`: The starting resid.
            - `resid2`: The ending resid.

            For example, a :class:`~pestifer.objs.deletion.Deletion` object can be used to define the range

        Yields
        ------
        :class:`~pestifer.molecule.residue.Residue`
            Yields residues within the specified range.
        """
        assert hasattr(rngrec,'chainID'), 'resrange requires a chainID'
        assert hasattr(rngrec,'resid1'), 'resrange requires a resid1'
        assert hasattr(rngrec,'resid2'), 'resrange requires a resid2'
        subR = self.get(lambda x: x.chainID == rngrec.chainID)
        subR.sort()
        r1 = rngrec.resid1
        r2 = rngrec.resid2
        R1 = subR.get(lambda x: x.resid == r1)
        if R1:
            R2 = subR.get(lambda x: x.resid == r2)
            if R2:
                idx1 = subR.index(R1)
                idx2 = subR.index(R2)
                assert idx2 >= idx1
                for j in range(idx1, idx2 + 1):
                    yield self[j]
        return []
    
    def do_deletions(self, Deletions: DeletionList):
        """
        Apply a list of deletions to the residue list.
        
        Parameters
        ----------
        Deletions : list of :class:`~pestifer.objs.deletion.Deletion`
            A list of deletion ranges to apply. Each deletion range should contain the chain ID, residue sequence numbers, and insertion codes for the start and end of the deletion range.
        """
        delete_us = []
        for d in Deletions:
            for dr in self.resrange(d):
                delete_us.append(dr)
        for d in delete_us:
            self.remove(d)

    def apply_segtypes(self):
        """
        Apply segment types to residues based on their residue names.
        This method uses the :attr:`~pestifer.core.labels.Labels.segtype_of_resname` mapping to assign segment types to
        residues based on their residue names. It updates the :attr:`~pestifer.molecule.residue.Residue.segtype` attribute of each residue
        in the residue list.
        """
        # logger.debug(f'residuelist:apply_segtypes {segtype_of_resname}')
        self.map_attr('segtype', 'resname', Labels.segtype_of_resname)

    def deletion(self, DL: DeletionList):
        """
        Remove residues from the residue list based on a :class:`~pestifer.objs.deletion.DeletionList`.
        This method iterates through the DeletionList, retrieves the residues to be deleted based on their
        chain ID, residue sequence numbers, and insertion codes, and removes them from the residue list
        
        Parameters
        ----------
        DL : :class:`~pestifer.objs.deletion.DeletionList`
            A list of deletion ranges to apply. Each deletion range should contain the chain ID, residue sequence numbers, and insertion codes for the start and end of the deletion range.
        """
        excised = []
        for d in DL:
            chain = self.filter(lambda x: x.chainID == d.chainID)
            r1 = chain.filter(lambda x: x.resid == d.resid1)[0]
            r2 = chain.filter(lambda x: x.resid == d.resid2)[0]
            for r in chain:
                if r1 <= r <= r2:
                    excised.append(r)
        for x in excised:
            self.remove(x)
            assert not x in self
        return excised

    def substitutions(self, SL:SubstitutionList):
        """
        Apply a list of substitutions to the residue list. This method iterates through the :class:`~pestifer.objs.substitution.SubstitutionList`, retrieves the residues to be substituted based on their chain ID, residue sequence numbers, and insertion codes, and replaces their residue names with the corresponding residue names from the substitution list. It also creates a new :class:`~pestifer.objs.seqadv.SeqadvList` for any resolved residues that are substituted.  Any residues that are deleted as a result of the substitutions are also returned in a list.

        Parameters
        ----------
        SL : :class:`~pestifer.objs.substitution.SubstitutionList`
            A list of substitutions to apply. Each substitution should contain the chain ID, residue sequence numbers,
            insertion codes, and the new residue sequence to substitute.
        
        Returns
        -------
        tuple : (:class:`~pestifer.objs.seqadv.SeqadvList`, list of :class:`~pestifer.molecule.residue.Residue`)
            A tuple containing a :class:`~pestifer.objs.seqadv.SeqadvList` of new sequence advancements for resolved residues and a list of residues that were deleted.
        """
        delete_us = []
        newseqadv = SeqadvList([])  # for holding single-residue changes for resolved residues
        for s in SL.data:
            subseq = s.subseq
            currsubidx = 0
            chain: ResidueList = self.get(lambda x: x.chainID == s.chainID)
            r1 = chain.filter(lambda x: x.resid == s.resid1)[0]
            r2 = chain.filter(lambda x: x.resid == s.resid2)[0]
            for r in chain.data:
                if r1 <= r <= r2:
                    if currsubidx < len(subseq):
                        resname = Labels.res_123[subseq[currsubidx].upper()]
                        if r.resolved:  # make a new seqadv for this mutation
                            newseqadv.append(Seqadv(
                                idCode='I doubt I ever use this',
                                typekey='user',
                                resname=r.resname,
                                chainID=r.chainID,
                                resid=ResID(r.resid),
                                dbRes=resname
                            ))
                        else:  # just change the residue name
                            r.resname = resname
                        currsubidx += 1
                    else:
                        delete_us.append(r)
            if currsubidx < len(subseq):
                # we have unsubstituted residue(s) left that must be inserted
                pass
        for r in delete_us:
            self.remove(r)
        return newseqadv, delete_us

    def cif_residue_map(self):
        """
        Create a mapping of residues by chain ID and residue sequence number.
        This method iterates through the residue list and creates a dictionary where the keys are chain IDs
        and the values are dictionaries mapping residue sequence numbers to Namespace objects containing
        the residue information.
        """
        result = {}
        for r in self.data:
            if hasattr(r, 'auth_asym_id'):  # this was constructed from atoms in a CIF file
                key = f'{r.chainID}:{r.resid.resid}'
                val = f'{r.auth_asym_id}:{r.auth_seq_id}{r.pdb_ins_code}'
                result[key] = val
        return result
    
    def apply_insertions(self, insertions: InsertionList):
        """
        Apply a list of insertions to the residue list.
        
        Parameters
        ----------
        insertions : list of :class:`~pestifer.objs.insertion.Insertion`
            A list of insertions to apply. Each insertion should contain the chain ID, residue sequence numbers,
            insertion codes, and the sequence of residues to insert.
        """
        for ins in insertions:
            chainID, resid = ins.chainID, ins.resid
            idx = self.iget(lambda x: x.chainID == chainID and x.resid == resid)
            segtype = self[idx].segtype
            logger.debug(f'insertion begins after {resid.resid} which is index {idx} in reslist, chain {chainID}')
            # add residues to residue list
            for olc in ins.sequence:
                resid.increment(by_seqnum=ins.integer_increment)
                shortcode = f'{chainID}:{olc}{resid.resid}'
                new_empty_residue = ResiduePlaceholder(shortcode)
                new_residue = Residue(new_empty_residue)
                new_residue.atoms = AtomList([])
                new_residue.segtype = segtype
                self.insert(idx, new_residue)
                logger.debug(f'insertion of new residue placeholder {shortcode}')

    def renumber(self, links: LinkList):
        """
        The possibility exists that empty residues added have resids that conflict with existing resids 
        on the same chain if those resids are in a different segtype (e.g., glycan).  This method will 
        privilege protein residues in such conflicts, and it will renumber non-protein residues, updating 
        any resid records in links to match the new resid.

        Parameters
        ----------
        links : :class:`~pestifer.objs.link.LinkList`
            A list of links that may contain residue sequence numbers that need to be updated.
        """
        protein_residues: ResidueList = self.get(lambda x: x.segtype == 'protein')
        if protein_residues is None or len(protein_residues) == 0: return
        min_protein_resid = min([x.resid for x in protein_residues.data])
        max_protein_resid = max([x.resid for x in protein_residues.data])
        non_protein_residues: ResidueList = ResidueList([])
        for p in self.data:
            if p.segtype != 'protein':
                non_protein_residues.append(p)
        logger.debug(f'There are {len(protein_residues)} (resids {min_protein_resid} to {max_protein_resid}) protein residues and {len(non_protein_residues)} non-protein residues')
        assert len(self) == (len(protein_residues) + len(non_protein_residues))
        non_protein_residues_in_conflict = ResidueList([])
        for np in non_protein_residues.data:
            tst = protein_residues.get(lambda x: x.chainID == np.chainID and x.resid == np.resid)
            if tst:
                non_protein_residues_in_conflict.append(np)
        for npc in non_protein_residues_in_conflict.data:
            non_protein_residues.remove(npc)
        logger.debug(f'There are {len(non_protein_residues_in_conflict)} non-protein residues with resid that conflict with protein residues')

        max_unused_resid = max([max([x.resid for x in protein_residues.data]), ResID(0) if len(non_protein_residues) == 0 else max([x.resid for x in non_protein_residues.data])]) + 1
        newtst = self.get(lambda x: x.resid == max_unused_resid)
        assert newtst in [None, []]  # None
        mapper_by_chain = {}
        for npc in non_protein_residues_in_conflict.data:
            c = npc.chainID
            if not c in mapper_by_chain:
                mapper_by_chain[c] = {}
            old_resid = npc.resid
            new = max_unused_resid
            mapper_by_chain[c][old_resid] = new
            max_unused_resid += 1
            npc.resid = new
            for a in npc.atoms.data:
                a.resid = new
            logger.debug(f'New resid: {c} {old_resid} -> {new}')
        if mapper_by_chain:
            logger.debug(f'Remapping resids in links')
            for l in links.data:
                if l.chainID1 in mapper_by_chain:
                    old_resid = l.resid1
                    if old_resid in mapper_by_chain[l.chainID1]:
                        l.resid1 = mapper_by_chain[l.chainID1][old_resid]
                        logger.debug(f' remapped left: {l.chainID1} {old_resid} -> {l.resid1}')
                if l.chainID2 in mapper_by_chain:
                    old_resid = l.resid2
                    if old_resid in mapper_by_chain[l.chainID2]:
                        l.resid2 = mapper_by_chain[l.chainID2][old_resid]
                        logger.debug(f' remapped right: {l.chainID2} {old_resid} -> {l.resid2}')

    def set_chainIDs(self, chainID):
        """
        Set the chain ID for all residues in the list to a specified value.
        
        Parameters
        ----------
        chainID : str
            The chain ID to set for all residues.
        """
        for r in self.data:
            r.set_chainID(chainID)

    def remap_chainIDs(self, the_map: dict):
        """
        Remap the chain IDs of residues in the list according to a provided mapping.

        Parameters
        ----------
        the_map : dict
            A dictionary mapping old chain IDs to new chain IDs.
        """
        for r in self.data:
            if r.chainID in the_map:
                r.set_chainID(the_map[r.chainID])

    def set_segname(self, segname: str):
        """
        Set the segment name for the residue.

        Parameters
        ----------
        segname : str
            The segment name to set.
        """
        for r in self.data:
            r.set_segname(segname)

    def state_bounds(self, state_func: Callable):
        """
        Get the state bounds for each residue in the list based on a state function.
        
        Parameters
        ----------
        state_func : function
            A function that takes a residue and returns its state bounds.
        
        Returns
        -------
        StateIntervalList
            A list of state intervals for each residue in the residue list.
        """

        return StateIntervalList.process_itemlist(self, state_func=state_func)

    def puniquify(self, attrs: list[str], stash_attr_name = 'ORIGINAL_ATTRIBUTES'):
        return super().puniquify(attrs, stash_attr_name)

Link.model_rebuild()
Patch.model_rebuild()
SSBond.model_rebuild()
Seqadv.model_rebuild()
StateInterval.model_rebuild()