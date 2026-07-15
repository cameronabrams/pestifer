#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Class for handling atoms.
"""
from __future__ import annotations

import logging

from functools import singledispatchmethod
from pidibble.pdbrecord import PDBRecord, PDBRecordDict
from pydantic import Field
from typing import ClassVar

from pestifer.util.stringthings import parse_filter_expression

from ..core.baseobj import BaseObj, BaseObjList
from ..objs.resid import ResID
from ..objs.ter import TerList
from ..psfutil.psfatom import PSFAtomList
from ..util.util import reduce_intlist

logger = logging.getLogger(__name__)

class Atom(BaseObj):
    """
    A class for handling atoms in molecular structures.
    This class represents an atom with various attributes such as serial number, name, residue name,
    chain ID, residue sequence number, insertion code, coordinates (x, y, z), occupancy, beta factor,
    element symbol, charge, and optional attributes like segment name, empty status, link status,
    record name, and author sequence ID, component ID, asym ID, and atom ID.
    """

    _required_fields = {'serial', 'name', 'altloc', 'resname', 'chainID', 'resid', 'x', 'y', 'z', 
                        'occ', 'beta', 'elem', 'charge'}
    """
    Required attributes for the Atom class.
    These attributes must be provided when creating an Atom instance.
    
    Attributes
    ----------
    serial : int
        Serial number of the atom.
    name : str 
        Name of the atom.
    altloc : str
        Alternate location identifier for the atom.
    resname : str
        Residue name to which the atom belongs.
    chainID : str
        Chain identifier for the atom.
    resid: ResID
        Residue ID of the atom.
    x : float
        X coordinate of the atom.
    y : float
        Y coordinate of the atom.
    z : float
        Z coordinate of the atom.
    occ : float
        Occupancy of the atom.
    beta : float
        Beta factor of the atom.
    elem : str
        Element symbol of the atom.
    charge : str
        Charge of the atom.
    """

    _optional_fields = {'segname', 'empty', 'link', 'recordname', 'auth_seq_id', 
                        'auth_comp_id', 'auth_asym_id', 'auth_atom_id', 'pdbx_pdb_ins_code',
                        'ORIGINAL_ATTRIBUTES'}
    """
    Optional attributes for the Atom class.
    These attributes can be provided when creating an Atom instance, but are not required.
    
    Attributes
    ----------
    segname : str
        Segment name to which the atom belongs. Defaults to the chain ID.
    empty : bool
        Indicates whether the atom is empty. Defaults to False.
    link : str
        Link status of the atom. Defaults to 'None'.
    recordname : str
        Record name for the atom. Defaults to 'ATOM'.
    auth_seq_id : int
        Author sequence ID for the atom.
    auth_comp_id : str
        Author component ID for the atom.
    auth_asym_id : str
        Author asym ID for the atom.
    auth_atom_id : str
        Author atom ID for the atom.
    """
    
    serial: int = Field(..., description="Serial number of the atom.")
    name: str = Field(..., description="Name of the atom.")
    altloc: str = Field(..., description="Alternate location identifier for the atom.")
    resname: str = Field(..., description="Residue name to which the atom belongs.")
    chainID: str = Field(..., description="Chain identifier for the atom.")
    resid: ResID = Field(..., description="Residue ID of the atom.")
    x: float = Field(..., description="X coordinate of the atom.")
    y: float = Field(..., description="Y coordinate of the atom.")
    z: float = Field(..., description="Z coordinate of the atom.")
    occ: float = Field(..., description="Occupancy of the atom.")
    beta: float = Field(..., description="Beta factor of the atom.")
    elem: str = Field(..., description="Element symbol of the atom.")
    charge: str = Field(..., description="Charge of the atom.")

    segname: str | None = Field(default=None, description="Segment name to which the atom belongs. Defaults to the chain ID.")
    empty: bool = Field(default=False, description="Indicates whether the atom is empty. Defaults to False.")
    link: str | None =  Field(default='None', description="Link status of the atom. Defaults to 'None'.")
    recordname: str | None = Field(default='ATOM', description="Record name for the atom. Defaults to 'ATOM'.")
    auth_seq_id: int | None = Field(default=None, description="Author sequence ID for the atom.")
    auth_comp_id: str | None = Field(default=None, description="Author component ID for the atom.")
    auth_asym_id: str | None = Field(default=None, description="Author asym ID for the atom.")
    auth_atom_id: str | None = Field(default=None, description="Author atom ID for the atom.")
    pdbx_pdb_ins_code: str | None = Field(default=None, description="PDB insertion code for the atom.")
    ORIGINAL_ATTRIBUTES: dict = Field(default_factory=dict, description="Dictionary to store original attributes of the atom instance.")

    _yaml_header: ClassVar[str] = 'atoms'
    """
    Header for YAML serialization of Atom objects.
    """

    _PDB_keyword: ClassVar[str] = 'ATOM'
    """
    Keyword used in PDB files to identify atom records.
    """

    _CIF_CategoryName: ClassVar[str] = 'atom_site'
    """
    Name used in mmCIF files to identify atom site records.
    """
    
    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        """
        Adapts the input to a dictionary format suitable for Atom instantiation.
        """
        if args and isinstance(args[0], PDBRecord):
            rec = args[0]
            # pidibble mmCIF records carry both a label (`residue`) and an author
            # (`residue_auth`) identity; PDB records carry only `residue` (author).
            is_cif = hasattr(rec, 'residue_auth')
            resname = rec.residue.resName
            occ = rec.occupancy
            beta = rec.tempFactor
            if resname == 'DUM':
                occ = 0.0
                beta = 0.0
            input_dict = {
                "serial": rec.serial,
                "name": rec.name,
                "resname": resname,
                "chainID": rec.residue.chainID,
                "x": rec.x,
                "y": rec.y,
                "z": rec.z,
                "occ": occ,
                "beta": beta,
                "elem": rec.element,
                # pidibble rectifies mmCIF pdbx_formal_charge to an int for charged
                # atoms; Atom.charge is a string (identity for the already-string PDB value).
                "charge": str(rec.charge),
                "altloc": rec.altLoc,
                "recordname": 'ATOM',
                "segname": rec.residue.chainID,
                "empty": False,
            }
            if is_cif:
                # mmCIF: label numbering is primary (chainID/resid), author numbering
                # is stashed as metadata (used later to map back to author chains and
                # to emit the cif_residue_map). pidibble rectifies purely-numeric chain
                # ids (e.g. '0') to ints, so coerce the chain-id fields back to strings.
                auth = rec.residue_auth
                input_dict['chainID'] = str(rec.residue.chainID)
                input_dict['segname'] = str(rec.residue.chainID)
                input_dict.update(
                    auth_seq_id=auth.seqNum,
                    auth_comp_id=auth.resName,
                    auth_asym_id=str(auth.chainID),
                    auth_atom_id=getattr(rec, 'auth_atom_id', None),
                    pdbx_pdb_ins_code=auth.iCode,
                )
                # non-polymer atoms have a blank label_seq_id; fall back to author number
                apparent_resseqnum = rec.residue.seqNum
                if apparent_resseqnum in ('', '.'):
                    apparent_resseqnum = auth.seqNum
                input_dict['resid'] = ResID(resseqnum=apparent_resseqnum, insertion=auth.iCode)
            else:
                input_dict['resid'] = ResID(resseqnum=rec.residue.seqNum, insertion=rec.residue.iCode)
            return input_dict
        return super()._adapt(*args, **kwargs)

    def shortcode(self) -> str:
        """
        Converts the Atom.Adapter object to a string representation.
        This method formats the attributes of the Atom.Adapter into a string.
        
        Returns
        -------
        str
            A string representation of the Atom.Adapter object.
        """
        return f"{self.serial}:{self.name}:{self.resname}:{self.chainID}:{self.resid.resid}:{self.x}:{self.y}:{self.z}:{self.occ}:{self.beta}:{self.elem}:{self.charge}"
    
    def pdb_line(self):
        """
        Returns a string representation of the atom in PDB format.
        This method formats the atom's attributes into a PDB line string.
        The line includes the record name, serial number, atom name, alternate location,
        residue name, chain ID, residue sequence number, insertion code, coordinates (x, y, z),
        occupancy, beta factor, element symbol, and charge.

        Returns
        -------
        str
            A formatted string representing the atom in PDB format.
        """
        pdbline='{:<6s}'.format(self.recordname)+\
                '{:5d}'.format(self.serial)+' '+\
                '{:<4s}'.format(' '+self.name if len(self.name)<4 else self.name)+\
                '{:1s}'.format(self.altloc)+\
                '{:<4s}'.format(self.resname)+\
                '{:1s}'.format(self.chainID)+\
                '{:>5s}'.format(self.resid.pdbresid)+'   '+\
                '{:8.3f}'.format(self.x)+\
                '{:8.3f}'.format(self.y)+\
                '{:8.3f}'.format(self.z)+\
                '{:6.2f}'.format(self.occ)+\
                '{:6.2f}'.format(self.beta)+\
                10*' '+'{:>2s}'.format(self.elem)+'{:2s}'.format(self.charge)
        return pdbline
    
    @singledispatchmethod
    def overwrite_position(self, *args):
        """
        Overwrites the position of this atom with the position of another atom.
        This method is a placeholder and should be overridden in subclasses.
        
        Parameters
        ----------
        *args : Any
            The arguments to overwrite the position. Should be an Atom object.
        
        Raises
        ------
        NotImplementedError
            If this method is called without being overridden in a subclass.
        """
        raise NotImplementedError("This method should be overridden in subclasses.")

    @overwrite_position.register(BaseObj)
    def _overwrite_position_from_Atom(self, other: BaseObj):
        """
        Overwrites the position of this atom with the position of another atom.
        This method updates the x, y, and z coordinates of this atom to match those of another atom.
        
        Parameters
        ----------
        other : Atom
            The atom whose position will be used to overwrite this atom's position.
        """
        if not isinstance(other, Atom):
            raise TypeError(f"Expected an Atom object, got {type(other)}")
        self.x=other.x
        self.y=other.y
        self.z=other.z

    @overwrite_position.register(dict)
    def _overwrite_position_from_dict(self, other: dict):
        """
        Overwrites the position of this atom with the position specified in a dictionary.
        The dictionary must contain keys 'x', 'y', and 'z' with the new coordinates.

        Parameters
        ----------
        other : dict
            A dictionary containing the new position for the atom.
        """
        self.x = other.get('x', self.x)
        self.y = other.get('y', self.y)
        self.z = other.get('z', self.z)

    @overwrite_position.register(float)
    def _overwrite_position_from_floats(self, x: float, y: float, z: float):
        """
        Overwrites the position of this atom with the position specified in a tuple.
        The tuple must contain three elements: (x, y, z) with the new coordinates.

        Parameters
        ----------
        x : float
            The new x coordinate.
        y : float
            The new y coordinate.
        z : float
            The new z coordinate.
        """
        self.x = x
        self.y = y
        self.z = z

class AtomList(BaseObjList[Atom]):
    """
    A class for handling lists of Atom objects.
    This class inherits from BaseObjList and provides methods to manage
    a list of Atom objects, including serialization, reserialization, and position overwriting.
    """

    def describe(self):
        return f'<AtomList with {len(self)} atoms>'

    @classmethod
    def from_pdb(cls, parsed: PDBRecordDict, model_id = None) -> "AtomList":
        """
        Create an AtomList from a PDBRecordDict.
        
        Parameters
        ----------
        parsed : PDBRecordDict
            The parsed PDB data containing atom records.
        
        Returns
        -------
        AtomList
            A new AtomList instance containing Atom objects created from the PDB data.
        """
        if Atom._PDB_keyword not in parsed:
            return cls([])
        return cls(
            [Atom(x) for x in parsed[Atom._PDB_keyword] if (model_id is None or x.model == model_id)]+
            [Hetatm(x) for x in parsed.get(Hetatm._PDB_keyword, []) if (model_id is None or x.model == model_id)]
        )

    @classmethod
    def from_cif(cls, parsed: PDBRecordDict) -> "AtomList":
        """
        Create an AtomList from a PDBRecordDict parsed from mmCIF.

        pidibble normalizes mmCIF `atom_site` into the same ATOM/HETATM records as PDB
        input (carrying both label and author residue identities), so this mirrors
        :meth:`from_pdb`. Rows are sorted by serial to recover atom_site file order,
        matching the previous mmcif-backed behavior.

        Parameters
        ----------
        parsed : PDBRecordDict
            The parsed mmCIF data.

        Returns
        -------
        AtomList
            A new AtomList instance containing Atom objects created from the mmCIF data.
        """
        atoms = cls(
            [Atom(x) for x in parsed.get(Atom._PDB_keyword, [])] +
            [Hetatm(x) for x in parsed.get(Hetatm._PDB_keyword, [])]
        )
        if len(atoms) == 0:
            return atoms
        atoms.sort(by=['serial'])
        return atoms

    def reserialize(self) -> AtomList:
        """
        Reserializes the AtomList by updating the serial numbers of each atom.
        This method assigns a new serial number to each atom in the list, starting from 1
        and incrementing for each atom. It also stores the original serial number in the
        `ORIGINAL_ATTRIBUTES` dictionary of each atom for reference.
        """
        serial = 1
        seen_per_chain: dict[str, set] = {}
        residshift_per_chain: dict[str, ResID] = {}
        current_resid_per_chain: dict[str, int | None] = {}
        for a in self.data:
            a.ORIGINAL_ATTRIBUTES['serial'] = a.serial
            a.ORIGINAL_ATTRIBUTES['resid'] = a.resid.copy(deep=True)
            a.serial = serial
            serial += 1
            chain = a.chainID
            if chain not in seen_per_chain:
                seen_per_chain[chain] = set()
                residshift_per_chain[chain] = ResID(0)
                current_resid_per_chain[chain] = None
            residshift = residshift_per_chain[chain]
            new_resid = a.resid + residshift
            # Only check collision on residue transitions, not between atoms of the same residue
            if new_resid.resseqnum != current_resid_per_chain[chain]:
                if new_resid.resseqnum in seen_per_chain[chain]:
                    residshift += ResID(9999)
                    new_resid = a.resid + residshift
                    residshift_per_chain[chain] = residshift
                    logger.debug(f'Atom {a.serial} chainID {chain} resid {a.resid} collides; shifting resseqnums by {residshift}')
                seen_per_chain[chain].add(new_resid.resseqnum)
                current_resid_per_chain[chain] = new_resid.resseqnum
            a.resid = new_resid
        return self
    
    def adjustSerials(self, Ters: TerList):
        """
        Adjusts the serial numbers of atoms in the AtomList based on the provided TerList.
        This method reduces the serial numbers of atoms in the AtomList by the number of
        ignored serials in the TerList. It updates the `ORIGINAL_` dictionary of each atom
        to store the original serial number before adjustment.

        Parameters
        ----------
        Ters : TerList
            A list of Ter objects containing serial numbers to be ignored (TER records in old-timey PDB files)
        """
        ignored_serials = [x.serial for x in Ters.data]
        if not ignored_serials:
            return
        logger.debug(f'These serials must be deleted: {ignored_serials}')
        ril = reduce_intlist([x.serial for x in self.data])
        logger.debug(f'Prior to ignore, serials populate {ril}')
        for a in self.data:
            try:
                n = next(x[0] for x in enumerate(ignored_serials) if x[1] > a.serial)
            except StopIteration:
                pass
            if n > 0:
                a.ORIGINAL_ATTRIBUTES['serial'] = a.serial
                a.serial -= n
                logger.debug(f'Atom orig serial {a.ORIGINAL_ATTRIBUTES["serial"]} to {a.serial}')

    def overwrite_positions(self, other: AtomList):
        """
        Overwrites the positions of atoms in this AtomList with the positions of atoms in another AtomList.
        This method iterates through both AtomLists, ensuring they are of equal length,
        and updates the position (x, y, z) of each atom in this AtomList to match the corresponding atom in the other AtomList.
        
        Parameters
        ----------
        other : AtomList
            The AtomList whose atom positions will be used to overwrite the positions in this AtomList.
        
        Raises
        ------
        AssertionError
            If the lengths of the two AtomLists are not equal, an assertion error is raised.
        """
        assert len(self) == len(other), 'Error: atom lists not equal length'
        for sa, oa in zip(self.data, other.data):
            sa.overwrite_position(oa)

    def apply_psf_attributes(self, psfatoms: PSFAtomList):
        """
        Applies attributes from a PSF atom list to the corresponding atoms in this AtomList.
        This method iterates through both AtomLists, ensuring they are of equal length,
        and updates the `resname` attribute of each atom in this AtomList to match the corresponding atom in the PSF atom list.
        
        Parameters
        ----------
        psfatoms : PSFAtomList
            The PSF atom list containing residue names to be applied.
        """
        for myatom, psfatom in zip(self.data, psfatoms.data):
            myatom.resname = psfatom.resname
            myatom.serial = psfatom.serial
            myatom.resid = psfatom.resid.copy(deep=True)
            myatom.segname = psfatom.segname

    def apply_inclusion_logics(self, inclusion_logics: list[str] = []) -> int:
        if len(inclusion_logics) == 0:
            return 0
        kept_atom_count = 0
        total_atom_count = len(self.data)
        keep_atoms = AtomList([])
        for expression in inclusion_logics:
            logger.debug(f'Applying atom inclusion logic: {expression}')
            filter_func = parse_filter_expression(expression)
            keep_atoms.extend(filter(filter_func, self.data))
        kept_atom_count = len(keep_atoms)
        logger.debug(f'Keeping {kept_atom_count} atoms.')
        if kept_atom_count > 0:
            self.data = keep_atoms
        return total_atom_count - kept_atom_count

    def apply_exclusion_logics(self, exclusion_logics: list[str] = []) -> int:
        if len(exclusion_logics) == 0:
            return 0
        seen_ids = set()
        all_ignored_atoms = AtomList([])
        for expression in exclusion_logics:
            logger.debug(f'Applying atom exclusion logic: {expression}')
            filter_func = parse_filter_expression(expression)
            for atom in filter(filter_func, self.data):
                if id(atom) not in seen_ids:
                    seen_ids.add(id(atom))
                    all_ignored_atoms.append(atom)
        logger.debug(f'Removing {len(all_ignored_atoms)} ignored atoms out of {len(self.data)} total atoms.')
        for atom in all_ignored_atoms:
            self.remove_instance(atom)
        return len(all_ignored_atoms)

class Hetatm(Atom):
    """
    A class for handling heteroatoms in molecular structures.
    This class inherits from the Atom class and represents heteroatoms with additional attributes.
    It includes the same attributes as Atom, but is specifically used for heteroatoms in PDB files.
    """

    _PDB_keyword: ClassVar[str] = 'HETATM'
    """
    Keyword used in PDB files to identify heteroatom records.
    """

    _yaml_header: ClassVar[str] = 'hetatoms'
    """
    Header for YAML serialization of Hetatm objects.
    """