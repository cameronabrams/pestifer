# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
This module defines classes and methods that facilitate reading and parsing of CHARMM RESI and MASS records from top, rtp, and str files from the force field. 
"""
from __future__ import annotations
from dataclasses import dataclass, field
import logging
import networkx as nx
import numpy as np

from collections import UserDict, UserList
from itertools import compress, batched

from ..scripters import PsfgenScripter
from ..util.stringthings import linesplit #, my_logger

logger=logging.getLogger(__name__)

@dataclass
class CharmmMass:
    """ 
    This class represents a single mass record from a CHARMM topology file.
    It contains the atom type, mass, and optionally the element symbol.
    The atom type is stored in uppercase, and the mass is stored as a float.
    """
    atom_type: str = ''
    """ The atom type, stored in uppercase. """
    atom_mass: float = 0.0
    """ The mass of the atom in AMU, stored as a float. """
    atom_element: str = '?'
    """ The element symbol, if provided, otherwise inferred from the atom type. """
    com: str = ''
    """ A comment associated with the mass record, if present. """

    @classmethod
    def from_card(cls, card: str) -> CharmmMass: 
        tok, com = linesplit(card)
        tokens = [t.strip() for t in tok.split()]
        atom_type = tokens[2].upper()
        atom_mass = float(tokens[3])
        atom_element = '?'
        try:
            atom_element = tokens[4]
        except:
            for c in atom_type[0]:
                if not c.isdigit():
                    atom_element = c.upper()
                    break
        return cls(
            atom_type=atom_type,
            atom_mass=atom_mass,
            atom_element=atom_element,
            com=com
        )

    def __str__(self):
        ans = f'{self.atom_type} {self.atom_mass:.4f}'
        if self.atom_element != '?':
            ans += f' {self.atom_element}'
        if self.com:
            ans += ' ' + self.com
        return ans

class CharmmMassList(UserList[CharmmMass]):

    def to_dict(self):
        return {m.atom_type: m for m in self}
    def from_dict(self, d: dict[str, CharmmMass]):
        self.clear()
        self.extend(d.values())

    @classmethod
    def from_cardlist(cls, cardlist: list[str]):
        mass_records = [CharmmMass.from_card(card) for card in cardlist]
        return cls(mass_records)

class CharmmMassDict(UserDict[CharmmMass]):
    """ 
    This class is a dictionary-like object that stores CharmmMassRecord objects.
    It allows for easy retrieval of mass records by atom type.
    
    Attributes
    ----------
    data : dict
        A dictionary mapping atom types to CharmmMassRecord objects.

    """
    def __init__(self, mList: list[CharmmMass]):
        self.data = {}
        for m in mList:
            self.data[m.atom_type] = m
    
    def to_list(self):
        return CharmmMassList(self.data.values())

    def from_list(self, other: CharmmMassList):
        self.data = {m.atom_type: m for m in other}

@dataclass
class CharmmAtom:
    """ 
    This class represents a single atom record from a CHARMM topology file.
    It contains the atom name, type, charge, mass, and optionally the element symbol.
    The atom name is stored in uppercase, and the type is also stored in uppercase.
    If the atom name starts with a digit, it is marked as an inpatch atom.
    """
    name: str = ''
    """ The name of the atom, stored in uppercase. """
    type: str = ''
    """ The type of the atom, stored in uppercase. """
    charge: float = 0.0
    """ The charge of the atom, stored as a float. """
    mass: float = 0.0
    """ The mass of the atom, stored as a float. """
    element: str = '?'
    """ The element symbol, if provided, otherwise set to '?'. """
    inpatch: bool = False
    """ A flag indicating whether the atom is part of a patch, set to True if the name starts with a digit. """

    @classmethod
    def from_card(cls, card: str) -> CharmmAtom:
        tokens = card.split()
        # logger.debug(f'Parsing CharmmAtom from card: [{card}]')
        return cls(
            name=tokens[1].upper(),
            type=tokens[2].upper(),
            charge=float(tokens[3]),
            inpatch=tokens[1][0].isdigit()
        )

    def __str__(self):
        ans = f'ATOM {self.name} {self.type} {self.charge:.4f}'
        if hasattr(self, 'mass') and self.mass > 0.0:
            ans += f' {self.mass:.4f}'
        if hasattr(self, 'element') and self.element != '?':
            ans += f' {self.element}'
        return ans

    def __eq__(self, other: CharmmAtom) -> bool:
        """ 
        Check if two CharmmAtom objects are equal. They are considered equal i
        f their names match, per the CHARMM specification.

        Parameters
        ----------
        other : CharmmAtom
            Another CharmmAtom object to compare against.

        Returns
        -------
        bool
            True if the names of the atoms match, False otherwise.
        """
        return self.name == other.name

    def set_mass(self, massdict: CharmmMassDict) -> CharmmAtom:
        """ 
        This method retrieves the mass record for the atom type from the provided CharmmMassDict object.
        If the atom type is found, it sets the mass and element
        attributes of the atom. If the atom type is not found, it logs an error and exits.
        
        Parameters
        ----------
        massdict : CharmmMassDict
            A CharmmMassDict object containing mass records for atom types.

        Raises
        ------
        ValueError
            If the atom type is not found in the CharmmMassDict object.
        TypeError
            If the masses parameter is not a CharmmMasses object.
        """

        if isinstance(massdict, CharmmMassDict):
            m: CharmmMass = massdict.get(self.type, None)
            if m is not None:
                self.mass = m.atom_mass
                self.element = m.atom_element
            else:
                raise ValueError(f'No mass record found for atom type {self.type}')
        else:
            raise TypeError(f'Cannot set mass of {self.name} ({self.type}) from {type(massdict)}')
        return self

class CharmmAtomList(UserList[CharmmAtom]):
    """ 
    This class is a list-like object that stores CharmmAtom objects.
    It allows for easy iteration over the atoms and provides methods to retrieve atoms by name, mass, serial number, and element.
    """
    def __init__(self, data: list[CharmmAtom]):
        self.data = data

    def __next__(self):
        """
        Iterates over the atoms in the list.
        """
        for i in range(len(self)):
            yield self[i]

    def to_dict(self, key: str = 'name'):
        return {getattr(atom, key): atom for atom in self.data}

    def from_dict(self, atom_dict: CharmmAtomDict):
        self.data = list(atom_dict.values())

    def get_atom(self, name: str) -> CharmmAtom | None:
        """ 
        This method searches for an atom by its name in the list of CharmmAtom objects.

        Parameters
        ----------
        name : str
            The name of the atom to retrieve.

        Returns
        -------
        CharmmAtom
            The CharmmAtom object with the specified name, or None if not found.
        """        
        L = [x.name for x in self.data]
        try:
            idx = L.index(name)
            return self.data[idx]
        except:
            return None

    def get_mass(self) -> float:
        """
        Return the sum of the masses of all atoms in the list.
        """
        return np.sum([a.mass for a in self.data])

    def get_mass_of_member(self, name: str) -> float:
        """ 
        This method searches for an atom by its name in the list of CharmmAtom objects
        and returns its mass. If the atom is not found or does not have a mass attribute, it returns 0.0.
        
        Parameters
        ----------
        name : str
            The name of the atom whose mass is to be retrieved.
        
        Returns
        -------
        float
            The mass of the atom in AMU, or 0.0 if the atom is not found or does not have a mass.
        """        
        a = self.get_atom(name)
        if a is not None:
            return a.mass
        return 0.0

    def get_serial(self, name: str) -> int:
        """ 
        Retrieve the serial number of an atom by its name.
        This method searches for an atom by its name in the list of CharmmAtom objects
        and returns its serial number. If the atom is not found or does not have a serial attribute, it returns -1.
        
        Parameters
        ----------
        name : str
            The name of the atom whose serial number is to be retrieved.
        
        Returns
        -------
        int
            The serial number of the atom, or -1 if the atom is not found or does not have a serial attribute.
        """
        a = self.get_atom(name)
        if a is not None:
            if hasattr(a, 'serial'):
                return getattr(a, 'serial')
        return -1

    def get_element(self, name: str) -> str:
        """ 
        Retrieve the element of an atom by its name.
        This method searches for an atom by its name in the list of CharmmAtom objects
        and returns its element symbol. If the atom is not found or does not have an element attribute, it returns '?'.
        
        Parameters
        ----------
        name : str
            The name of the atom whose element is to be retrieved.
        
        Returns
        -------
        str
            The element symbol of the atom, or '?' if the atom is not found or does not have an element attribute.
        """
        a = self.get_atom(name)
        if a is not None:
            if hasattr(a, 'element'):
                if a.element == '?':
                    raise ValueError(f'{a.name} has no element?')
                # logger.debug(f'returning {a.element} for {name}')
                return a.element
        else:
            logger.debug(f'Could not find atom with name {name}')
        logger.error(f'{name} has no element?')
        return '?'

    def append(self, atom: CharmmAtom):
        """ 
        This method appends a CharmmAtom object to the caller.
        If the atom does not have a serial number, it assigns a new serial number based on the current length of the list.
        
        Parameters
        ----------
        atom : CharmmAtom
            The CharmmAtom object to append to the list.
        """
        if not hasattr(atom, 'serial') or getattr(atom, 'serial') is None:
            setattr(atom, 'serial', len(self) + 1)
        super().append(atom)

    def set_masses(self, massdict: CharmmMassDict):
        """ 
        This method iterates through all CharmmTopAtom objects in the list and sets their masses
        using the provided CharmmMassDict object. It calls the `set_mass` method of each atom.

        Parameters
        ----------
        massdict : CharmmMassDict
            A CharmmMassDict object containing mass records for atom types.
        """        
        for a in self.data:
            # logger.debug(f'setting mass for {type(a)} \'{str(a)}\'')
            a.set_mass(massdict)

class CharmmAtomDict(UserDict[str, CharmmAtom]):
    """ 
    This class represents a dictionary mapping atom names to their corresponding CharmmAtom objects.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __setitem__(self, key: str, value: CharmmAtom):
        if not isinstance(value, CharmmAtom):
            raise TypeError(f"Expected CharmmAtom, got {type(value)}")
        super().__setitem__(key, value)

    def __getitem__(self, key: str) -> CharmmAtom:
        return super().__getitem__(key)

    def to_list(self) -> CharmmAtomList:
        return CharmmAtomList(self.values())
    
    def from_list(self, atom_list: CharmmAtomList):
        self.clear()
        for atom in atom_list:
            self[atom.name] = atom

@dataclass
class CharmmBond:
    """ 
    This class represents a single bond record from a CHARMM topology file.
    It contains the names of the two atoms that form the bond and the degree of the bond (single, double, or triple).
    The names are stored in uppercase, and the degree is stored as an integer.
    Because multiple bonds are stored on a line, input is managed by CharmmBondList.
    """

    name1: str = ''
    """ The name of the first atom in the bond, stored in uppercase. """
    name2: str = ''
    """ The name of the second atom in the bond, stored in uppercase. """
    degree: int = 1
    """ The degree of the bond, which can be 1 (single), 2 (double), or 3 (triple). Defaults to 1. """
    length: float = 0.0
    """ Length of the bond in Angstroms; only relevant within an IC """

    def __eq__(self, other: CharmmBond) -> bool:
        """ 
        This method compares the names of the atoms in the bond, allowing for the order of the names to be reversed.
        It returns True if the names match in either order, and False otherwise.
        
        Parameters
        ----------
        other : CharmmBond
            Another CharmmBond object to compare against.
        
        Returns
        -------
        bool
            True if the names of the atoms in the bond match in either order, False otherwise.
        """
        c1 = self.name1 == other.name1 and self.name2 == other.name2
        c2 = self.name1 == other.name2 and self.name2 == other.name1
        return c1 or c2
    
class CharmmBondList(UserList[CharmmBond]):
    """ 
    This class is a list-like object that stores CharmmBond objects.
    It allows for easy iteration over the bonds and provides methods to retrieve bonds by atom names.
    
    Attributes
    ----------
    data : list
        A list of CharmmBond objects representing the bonds in a CHARMM topology file.
    """
    def __init__(self, data: list[CharmmBond]):
        self.data = data

    @classmethod
    def from_card(cls, card: str, degree: int = 1) -> "CharmmBondList":
        """Create a CharmmBondList from a CHARMM bond record."""
        result = CharmmBondList([])
        toks = [x.strip() for x in card.split()][1:]
        for p in batched(toks, 2):
            result.append(CharmmBond(*p, degree=degree))
        return result

@dataclass
class CharmmAngle:
    """ 
    This class represents a single angle record from a CHARMM topology file.
    It contains the names of the three atoms that form the angle and the angle in degrees.
    The names are stored in uppercase, and the angle is stored as a float.
    """
    name1: str = ''
    """ The name of the first atom in the angle, stored in uppercase. """
    name2: str = ''
    """ The name of the second atom in the angle, stored in uppercase. """
    name3: str = ''
    """ The name of the third atom in the angle, stored in uppercase. """
    degrees: float = 0.0
    """ The angle in degrees; only relevant in an IC """

    def __eq__(self, other: CharmmAngle) -> bool:
        """ 
        This method compares the names of the atoms in the angle, allowing for the order of the names to be reversed.
        It returns True if the names match in either order, and False otherwise.
        
        Parameters
        ----------
        other : CharmmAngle
            Another CharmmAngle object to compare against.
        
        Returns
        -------
        bool
            True if the names of the atoms in the angle match in either order, False otherwise.
        """
        r1 = self.name2 == other.name2
        c1 = self.name1 == other.name1 and self.name3 == other.name3
        c2 = self.name1 == other.name3 and self.name3 == other.name1
        return r1 and (c1 or c2)

class CharmmAngleList(UserList[CharmmAngle]):
    """ 
    This class is a list-like object that stores CharmmAngle objects.
    It allows for easy iteration over the angles and provides methods to retrieve angles by atom names.

    Attributes
    ----------
    data : list
        A list of CharmmAngle objects representing the angles in a CHARMM topology file.
    """
    def __init__(self, data: list[CharmmAngle]):
        self.data = data

    @classmethod
    def from_card(cls, card: str) -> "CharmmAngleList":
        result = cls([])
        toks = [x.strip() for x in card.split()][1:]
        for p in batched(toks, 3):
            result.append(CharmmAngle(*p))
        return result

@dataclass
class CharmmDihedral:
    """ 
    This class represents a single dihedral record from a CHARMM topology file.
    It contains the names of the four atoms that form the dihedral and the dihedral type (proper or improper).
    The names are stored in uppercase, and the dihedral type is stored as a string.
    """

    name1: str = ''
    """ The name of the first atom in the dihedral, stored in uppercase. """
    name2: str = ''
    """ The name of the second atom in the dihedral, stored in uppercase. """
    name3: str = ''
    """ The name of the third atom in the dihedral, stored in uppercase. """
    name4: str = ''
    """ The name of the fourth atom in the dihedral, stored in uppercase. """
    dihedral_type: str = 'proper'
    """ The type of the dihedral, which can be 'proper' or 'improper'. Defaults to 'proper'. """
    degrees: float = 0.0
    """ The angle in degrees; only relevant in an IC. """

    def __eq__(self, other: CharmmDihedral) -> bool:
        """ 
        This method compares the names of the atoms in the dihedral and the dihedral type.
        It returns True if the names match in the order they are defined and the dihedral types are the same, and False otherwise.
        
        Parameters
        ----------
        other : CharmmDihedral
            Another CharmmDihedral object to compare against.
        
        Returns
        -------
        bool
            True if the names of the atoms in the dihedral match in the order they are defined and the dihedral types are the same, False otherwise.
        """
        l1 = [self.name1, self.name2, self.name3, self.name4]
        l2 = [other.name1, other.name2, other.name3, other.name4]
        if self.dihedral_type != other.dihedral_type:
            return False
        return l1 == l2

class CharmmDihedralList(UserList[CharmmDihedral]):
    """ 
    This class is a list-like object that stores CharmmDihedral objects.
    It allows for easy iteration over the dihedrals and provides methods to retrieve dihedrals by atom names.

    Attributes
    ----------
    data : list
        A list of CharmmDihedral objects representing the dihedrals in a CHARMM topology file.
    """
    def __init__(self, data: list[CharmmDihedral]):
        self.data = data

    @classmethod
    def from_card(cls, card: str, dihedral_type: str = 'proper') -> "CharmmDihedralList":
        dihedrals = CharmmDihedralList([])
        toks = [x.strip() for x in card.split()][1:]
        for p in batched(toks, 4):
            dihedrals.append(CharmmDihedral(*p, dihedral_type=dihedral_type))
        return dihedrals

@dataclass
class CharmmIC:
    """ 
    This class represents a single internal coordinate (IC) record from a CHARMM topology file.
    It contains the names of the four atoms that define the internal coordinate, the lengths of the bonds, the angles, and the dihedral angle.
    The names are stored in uppercase, and the lengths, angles, and dihedral angle are stored as floats.
    """

    card: str = ''
    """ The unprocessed string for the IC """
    atoms: CharmmAtomList = field(default_factory=CharmmAtomList)
    """ A list of the names of the four atoms that define the internal coordinate, stored in uppercase. """
    bonds: CharmmBondList = field(default_factory=CharmmBondList)
    """ A CharmmBondList containing the bonds defined by the internal coordinate. """
    angles: CharmmAngleList = field(default_factory=CharmmAngleList)
    """ A CharmmAngleList containing the angles defined by the internal coordinate. """
    dihedral: CharmmDihedral = field(default_factory=CharmmDihedral)
    """ A CharmmDihedral object representing the dihedral angle defined by the internal coordinate. """
    dihedral_type: str = 'proper'
    """ The type of the dihedral, which can be 'proper' or 'improper'. Defaults to 'proper'. """
    empty: bool = False
    """ A flag indicating whether the internal coordinate is empty (i.e., missing data). Defaults to False. """

    @classmethod
    def from_card(cls, card: str) -> "CharmmIC":
        empty = False
        ICstring, dummy = linesplit(card)
        toks = [x.strip() for x in ICstring.split()]
        data = np.array(list(map(float, toks[5:])))
        if data[0] == 0.0 or data[1] == 0.0 or data[3] == 0.0 or data[4] == 0.0:
            # logger.debug(f'{toks[1:5]}: missing ic data: {data}')
            empty = True
        atoms = toks[1:5]
        if atoms[2].startswith('*'):
            atoms[2] = atoms[2][1:]
            dihedral_type = 'improper'
        else:
            dihedral_type = 'proper'
        dihedral = CharmmDihedral(name1=atoms[0], name2=atoms[1], name3=atoms[2], name4=atoms[3], dihedral_type=dihedral_type, degrees=float(toks[7]))
        if dihedral_type == 'proper':
            b0: CharmmBond = CharmmBond(name1=atoms[0], name2=atoms[1], length=float(toks[5]))
            b1: CharmmBond = CharmmBond(name1=atoms[2], name2=atoms[3], length=float(toks[9]))
            bonds = CharmmBondList([b0, b1])
            a1: CharmmAngle = CharmmAngle(name1=atoms[0], name2=atoms[1], name3=atoms[2], degrees=float(toks[6]))
            a2: CharmmAngle = CharmmAngle(name1=atoms[1], name2=atoms[2], name3=atoms[3], degrees=float(toks[8]))
            angles = CharmmAngleList([a1, a2])
        else:
            b0 = CharmmBond(name1=atoms[0], name2=atoms[2], length=float(toks[5]))
            b1 = CharmmBond(name1=atoms[2], name2=atoms[3], length=float(toks[9]))
            bonds = CharmmBondList([b0, b1])
            a1 = CharmmAngle(name1=atoms[0], name2=atoms[2], name3=atoms[3], degrees=float(toks[6]))
            a2 = CharmmAngle(name1=atoms[1], name2=atoms[2], name3=atoms[3], degrees=float(toks[8]))
            angles = CharmmAngleList([a1, a2])
        return cls(atoms=atoms, bonds=bonds, angles=angles, dihedral=dihedral, dihedral_type=dihedral_type, empty=empty)

    def __str__(self):
        ans = ','.join([f'{i+1}{self.atoms[i]}' for i in range(len(self.atoms))])
        ans += ':'
        ans += ','.join([f'{b.length:.4f}' for b in self.bonds])
        ans += ':'
        ans += ','.join([f'{a.degrees:.4f}' for a in self.angles])
        ans += ':'
        ans += f'{self.dihedral.degrees}-{self.dihedral.dihedral_type}'
        if self.empty:
            ans += '-empty'
        return ans

class CharmmICList(UserList[CharmmIC]):
    """
    This class represents a list of internal coordinate records from a CHARMM topology file.
    """
    pass

@dataclass
class CharmmDelete:
    """ 
    This class represents a single delete atom record from a CHARMM topology file.
    It contains the atom name to be deleted and an optional comment.
    The atom name is stored in uppercase, and the comment is stored as a string.
    """    

    atomname: str = ''
    """ The name of the atom to be deleted, stored in uppercase. """
    comment: str = ''
    """ An optional comment associated with the delete atom record. """
    raw_string: str = ''
    """ The raw string representation of the delete atom record. """
    error_code: int = 0
    """ An error code indicating the status of the delete atom record. Defaults to 0. """

    @classmethod
    def from_card(cls, card: str) -> 'CharmmDelete':
        instance = cls()
        instance.raw_string = card
        instance.error_code = 0
        card, instance.comment = linesplit(card)
        toks = [x.strip() for x in card.split()]
        if len(toks) < 3 or not toks[0].upper().startswith('DELETE'):
            instance.error_code = -1
        instance.atomname = toks[2].upper()
        return instance

    def __str__(self):
        return f'DELETE ATOM {self.atomname} ! {self.comment}'

class CharmmDeleteList(UserList[CharmmDelete]):
    """
    This class represents a list of delete atom records from a CHARMM topology file.
    """
    pass

@dataclass
class CharmmResi:
    """ 
    This class represents a single residue record (RESI) from a CHARMM topology file.
    It contains the residue name, charge, synonym, and a list of atoms, bonds,
    internal coordinates (ICs), and delete atoms associated with the residue.
    The residue name is stored in uppercase, and the charge is stored as a float.
    """

    key: str = ''
    """ The key for the residue record, typically 'RESI' or 'PRES'. """
    blockstring: str = ''
    """ The raw string representation of the residue record. """
    resname: str = ''
    """ The name of the residue, stored in uppercase. """
    charge: float = 0.0
    """ The charge of the residue, stored as a float. """
    mass: float = 0.0
    """ The mass of the residue, stored as a float; computed after all mass records processed """
    synonym: str = ''
    """ A synonym or comment associated with the residue. """
    metadata: dict = field(default_factory=dict)
    """ A dictionary containing metadata associated with the residue. """
    atoms: CharmmAtomList = field(default_factory=CharmmAtomList)
    """ A list of CharmmAtom objects representing the atoms in the residue. """
    atoms_in_group: dict[str, CharmmAtomList] = field(default_factory=dict)
    """ A dictionary mapping atom group indices to lists of CharmmAtom objects. """
    bonds: CharmmBondList = field(default_factory=CharmmBondList)
    """ A list of CharmmBond objects representing the bonds in the residue. """
    ICs: CharmmICList = field(default_factory=CharmmICList)
    """ A list of CharmmTopIC objects representing the internal coordinates in the residue. """
    Delete: CharmmDeleteList = field(default_factory=CharmmDeleteList)
    """ A list of CharmmTopDelete objects representing the delete atom records in the residue. """
    error_code: int = 0
    """ An error code indicating the status of the residue record. Defaults to 0. """
    ispatch: bool = False
    """ A flag indicating whether the residue contains atoms that are part of a patch. Defaults to False. """
    source_file: str = ''
    """ The source file from which the residue was constructed. """

    @staticmethod
    def name_from_blockstring(blockstring: str):
        lines = [x.strip() for x in blockstring.split('\n')]
        titledata, _ = linesplit(lines[0])
        tctokens = titledata.split()
        return tctokens[1] if len(tctokens) > 1 else ''

    @classmethod
    def from_blockstring(cls, blockstring: str,  metadata: dict = {}) -> 'CharmmResi':
        # logger.debug(f'Parsing CharmmResi from blockstring: {blockstring[:30]}...')

        error_code = 0
        lines = [x.strip() for x in blockstring.split('\n')]
        titledata, titlecomment = linesplit(lines[0])
        tctokens = titledata.split()
        key = tctokens[0].upper()
        resname = tctokens[1]
        charge = float(tctokens[2])
        synonym = titlecomment.strip()
        metadata = metadata

        didx = 1
        while lines[didx].startswith('!'): didx += 1
        cards = [x.upper() for x in lines[didx:]]
        # logger.debug(f'first post-title card {cards[0]}')
        comments = []
        datacards = []
        for c in cards:
            rawcarddata, cardcomment = linesplit(c)
            carddata = rawcarddata.lstrip().upper()  # no indent, no case-sensitivity
            if len(carddata) > 0: datacards.append(carddata)
            if len(cardcomment) > 0: comments.append(cardcomment)
        # logger.debug(f'{self.resname}: {len(datacards)} datacards and {len(comments)} comments')
        isatomgroup = [d.startswith('ATOM') or d.startswith('GROU') for d in datacards]
        atomgroupcards = list(compress(datacards, isatomgroup))
        # logger.debug(f'{self.resname}: {len(atomgroupcards)} atom group cards')
        # logger.debug(f'{atomgroupcards[0:5]}')
        atoms = CharmmAtomList([])
        atoms_in_group: dict[int, CharmmAtomList] = {}
        g = 0
        for card in atomgroupcards:
            # logger.debug(f'{card}')
            if card.startswith('ATOM'):
                a = CharmmAtom.from_card(card)
                atoms.append(a)
                if not g in atoms_in_group:
                    atoms_in_group[g] = CharmmAtomList([])
                atoms_in_group[g].append(a)
            elif card.startswith('GROU'):
                g += 1
        ispatch = False
        for a in atoms:
            if hasattr(a, 'inpatch'):
                if a.inpatch:
                    ispatch = True
                    break
        # logger.debug(f'{len(self.atoms)} atoms processed in {len(self.atoms_in_group)} groups')
        atomdict = atoms.to_dict(key='name')
        isbond = [d.startswith('BOND') for d in datacards]
        bondcards = compress(datacards, isbond)

        bonds = CharmmBondList([])
        for card in bondcards:
            bonds.extend(CharmmBondList.from_card(card))
        bondcards = None
        isdouble = [d.startswith('DOUB') for d in datacards]
        doublecards = compress(datacards, isdouble)
        for card in doublecards:
            bonds.extend(CharmmBondList.from_card(card, degree=2))
        istriple = [d.startswith('TRIP') for d in datacards]
        triplecards = compress(datacards, istriple) 
        for card in triplecards:
            bonds.extend(CharmmBondList.from_card(card, degree=3))
        isIC = [d.startswith('IC') for d in datacards]
        ICcards = compress(datacards, isIC)
        ICs = CharmmICList([])
        ic_atom_names = []
        for card in ICcards:
            IC = CharmmIC.from_card(card)
            if any([x not in atomdict.keys() for x in IC.atoms]):
                error_code = -5
            else:
                ic_atom_names.extend(IC.atoms)
                ICs.append(IC)
        isDelete = [d.startswith('DELETE') for d in datacards]
        deletecards = compress(datacards, isDelete)
        Delete = CharmmDeleteList([])
        for card in deletecards:
            D = CharmmDelete.from_card(card)
            Delete.append(D)

        return cls(key=key, blockstring=blockstring, resname=resname, charge=charge, synonym=synonym, metadata=metadata, atoms=atoms,atoms_in_group=atoms_in_group, bonds=bonds, ICs=ICs, Delete=Delete, error_code=error_code, ispatch=ispatch)

    def copy_ICs_from(self, other: CharmmResi):
        """ 
        This method copies the internal coordinates from another CharmmResi object to the current object.
        It checks if the atoms in the ICs exist in the current object's atom dictionary.
        If any atom in the IC is not found in the current object's atom dictionary, it logs an error and sets the error code to -5.
        
        Parameters
        ----------
        other : CharmmResi
            Another CharmmResi object from which to copy the internal coordinates.
        """
        self.IC = []
        for ic in other.IC.data:
            ic_card = ic.card
            IC = CharmmIC(ic_card)
            if any([x not in self.atomdict.keys() for x in IC.atoms]):
                logger.error(f'IC {ic_card} has atoms not in {self.resname}: {IC.atoms}')
                self.error_code = -5
            else:
                self.IC.append(IC)

    def set_mass(self, massdict: CharmmMassDict):
        """ 
        This method iterates through all CharmmTopAtom objects in the residue and sets their masses
        using the provided CharmmMasses object. It calls the :meth:`~pestifer.charmmff.charmmfftop.CharmmTopAtom.set_mass` method of each atom.

        Parameters
        ----------
        masses : CharmmMasses
            A CharmmMasses object containing mass records for atom types.
        """
        self.atoms.set_masses(massdict)
        self.mass = self.atoms.get_mass()

    def num_atoms(self):
        """ 
        Return the number of atoms in the residue.
        This method returns the length of the atoms list in the residue.
        
        Returns
        -------
        int
            The number of atoms in the residue.
        """
        return len(self.atoms)

    def formula(self):
        """ 
        This method constructs an empirical chemical formula string based on the elements present in the atoms of the residue.
        It counts the occurrences of each element and formats the string according to the standard chemical notation.
        The elements are sorted in a predefined order (C, H, N, O, P, and then others).

        Returns
        -------
        str
            The empirical chemical formula of the residue, formatted as a string. If no elements are found, it returns an empty string.
        """
        fdict = {}
        sortorder = 'CHNOP'
        for a in self.atoms:
            if a.element != '?':
                if not a.element in fdict:
                    fdict[a.element] = 0
                fdict[a.element] += 1
        if fdict:
            retstr = ''
            for e in sortorder:
                if e in fdict:
                    n = '' if fdict[e] == 1 else str(fdict[e])
                    retstr += f'{e}{n}'
            for k, v in fdict.items():
                if not k in sortorder:
                    n = '' if v == 1 else str(v)
                    retstr += f'{k}{n}'
        return retstr

    def to_graph(self, includeH: bool = True):
        """ 
        This method creates a :class:`networkx.Graph` instance from the bonds in the residue.
        It adds edges between atoms based on the bonds, and assigns the element of each atom as a node attribute.
        If parameter `includeH` is True, hydrogen atoms are included in the graph; otherwise, they are excluded.
        
        Parameters
        ----------
        includeH : bool, optional
            A flag indicating whether to include hydrogen atoms in the graph. Defaults to True.
        
        Returns
        -------
        networkx.Graph
            A networkx graph representing the residue, with atoms as nodes and bonds as edges. Each node has an 'element' attribute representing the element of the atom.
        """
        g = nx.Graph()
        for b in self.bonds.data:
            node1 = b.name1
            node2 = b.name2
            # logger.debug(f'bond {node1} {node2}')
            element1 = self.atoms.get_element(node1)
            element2 = self.atoms.get_element(node2)
            assert element1 != '?', f'no element for atom {node1} in {self.resname}'
            assert element2 != '?', f'no element for atom {node2} in {self.resname}'
            if element1 == 'H' or element2 == 'H':
                if includeH:
                    g.add_edge(node1, node2)
                    g.nodes[node1]['element'] = element1
                    g.nodes[node2]['element'] = element2
            else:
                g.add_edge(node1, node2)
                g.nodes[node1]['element'] = element1
                g.nodes[node2]['element'] = element2
        return g
    
    def lipid_annotate(self):
        """ 
        This method performs lipid-specific annotation based on the metadata of the residue.
        It checks the `streamID` and `substreamID` in the metadata to determine
        the type of lipid and calls the appropriate annotation method.
        The available annotation methods are:
        - `sterol_annotate`: for cholesterol lipids
        - `detergent_annotate`: for detergent lipids
        - `model_annotate`: for model lipids
        - `generic_lipid_annotate`: for generic lipids
        The method initializes the `annotation` attribute with heads, tails, and shortest paths based on the lipid type.
        """
        self.annotation = {}
        m = self.metadata
        logger.debug(f'metadata {m}')
        if m['streamID'] == 'lipid' and m['substreamID'] == 'cholesterol':
            self.sterol_annotate()
        elif m['streamID'] == 'lipid' and m['substreamID'] == 'detergent':
            self.detergent_annotate()
        elif m['streamID'] == 'lipid' and m['substreamID'] == 'model':
            self.model_annotate()
        else:
            self.generic_lipid_annotate()
    
    def model_annotate(self):
        """ 
        This method initializes the `annotation` attribute with empty heads, tails, and shortest paths.
        It is a placeholder and does not perform any specific operations.
        """
        self.annotation = {}
        self.annotation['heads'] = []
        self.annotation['tails'] = []
        self.annotation['shortest_paths'] = {}

    def sterol_annotate(self):
        """ 
        This method identifies the head and tail of the sterol based on the carbon atoms in the residue.
        It constructs a graph representation of the residue and finds the carbon atom to which the hydroxyl oxygen is bound as the head.
        The tail is determined as the carbon atom furthest away from the head.
        The method initializes the `annotation` attribute with heads, tails, and shortest paths based on the sterol structure.
        """
        self.annotation = {}
        G = self.to_graph(includeH=False)
        heads = []
        for a in G:
            # head atom is the carbon to which the hydroxyl O is bound
            for n in nx.neighbors(G, a):
                if G.nodes[n]['element'] == 'O':
                    heads = [a]
                    break
        paths = []
        for a in G.__iter__():
            paths.append(len(nx.shortest_path(G, a, heads[0])))
        # tail is atom furthest away from head
        paths = np.array(paths)
        l_idx = np.argsort(paths)[-1]
        tails = [list(G.__iter__())[l_idx]]
        logger.debug(f'heads {heads} tails {tails}')
        self.annotation['heads'] = heads
        self.annotation['tails'] = tails

    def detergent_annotate(self):
        """ 
        This method identifies the head and tail of the detergent based on the carbon atoms in the residue.
        It constructs a graph representation of the residue and finds the carbon atom to which the hydroxyl oxygen is bound as the head.
        The tail is determined as the carbon atom furthest away from the head.
        The method initializes the `annotation` attribute with heads, tails, and shortest paths based on the detergent structure.
        """
        self.annotation = {}
        G = self.to_graph(includeH=False)
        # find the bond which, if deleted, produces two fragments such that
        # 1. one fragment contains only carbons (tail)
        # 2. the other fragment is as small as possible and contains at least one O (head)
        minheadsize = 1000
        mindata = {}
        for f in nx.bridges(G):
            logger.debug(f'f {f}')
            g = G.copy()
            g.remove_edge(*f)
            c = list(nx.connected_components(g))
            # if len(c)!=2: continue
            e = [[G.nodes[x]['element'] for x in y] for y in c]
            logger.debug(f'e {e}')
            allc = [all([x == 'C' for x in y]) for y in e]
            haso = [any([x == 'O' for x in y]) for y in e]
            logger.debug(f'allc {allc}')
            logger.debug(f'haso {haso}')
            if any(allc) and any(haso):
                tailidx = allc.index(True)
                headidx = 1 - tailidx
                logger.debug(f'tailidx {tailidx} headidx {headidx}')
                tailsize = len(e[tailidx])
                headsize = len(e[headidx])
                logger.debug(f'tailsize {tailsize} headsize {headsize}')
                if headsize < minheadsize:
                    minheadsize = headsize
                    mindata = dict(headg=list(c)[headidx], tailg=list(c)[tailidx])
        logger.debug(f'mindata {mindata}')
        headg = list(mindata['headg'])
        tailg = list(mindata['tailg'])
        # head reference is heaviest atom in head
        headmasses = np.array([self.atoms.get_mass_of_member(n) for n in headg])
        headatom = headg[np.argmax(headmasses)]
        logger.debug(f'headmasses {headmasses} headatom {headatom}')

        self.annotation['heads'] = [headatom]

        maxpathlength = -1
        tailatom = None
        for ta in tailg:
            path = nx.shortest_path(G, source=headatom, target=ta)
            pathlength = len(path)
            if pathlength > maxpathlength:
                maxpathlength = pathlength
                tailatom = ta
        assert tailatom != None
        self.annotation['tails'] = [tailatom]
        WG = G.copy()
        self.annotation['shortest_paths'] = {
            head: {
                tail: len(nx.shortest_path(WG, source=head, target=tail)) for tail in self.annotation['tails']
            } for head in self.annotation['heads']
        }

        logger.debug(f'annotation {self.annotation}')

    def generic_lipid_annotate(self):
        """ 
        This method identifies the heads and tails of the lipid based on the carbon atoms in the residue.
        It constructs a graph representation of the residue and finds carbon chains (tails) by removing edges
        that, when removed, produce two fragments: one containing only carbons (tail) and the other containing at least one oxygen (head).
        The method initializes the `annotation` attribute with heads, tails, and shortest paths based on the lipid structure.
        It handles cases where there are two tails (as in standard lipids) or multiple chains (as in complex lipids).
        """
        self.annotation = {}
        G = self.to_graph(includeH=False)
        WG = G.copy()

        # find all carbon chains
        ctails = []
        cbonds = []
        hastails = True
        while hastails:
            hastails = False
            maxtailsize = -1
            result = {}
            for f in nx.bridges(WG):
                g = WG.copy()
                g.remove_edge(*f)
                S = [g.subgraph(c).copy() for c in nx.connected_components(g)]
                c = list(nx.connected_components(g))
                if len(c) != 2: continue
                e = [[G.nodes[x]['element'] for x in y] for y in c]
                # logger.debug(f'e {e}')
                allc = [all([x == 'C' for x in y]) for y in e]
                # haso=[any([x=='O' for x in y]) for y in e]
                # logger.debug(f'allc {allc}')
                # logger.debug(f'haso {haso}')
                if any(allc):
                    tail = c[allc.index(True)]
                    tailsize = len(tail)
                    # logger.debug(f'tail {tail} tailsize {tailsize}')
                    if tailsize > maxtailsize:
                        result = dict(tail=tail, nontailg=S[allc.index(False)], bond=f)
                        tailatoms = list(tail)
                        maxtailsize = tailsize
            if result:
                hastails = True
                # ignore methyls
                if len(tailatoms) > 1:
                    ctails.append(tailatoms)
                    cbonds.append(result['bond'])
                    logger.debug(f'tailatoms {tailatoms}')
                WG = result['nontailg']

        logger.debug(f'ctails {ctails}')
        logger.debug(f'cbonds {cbonds}')
        WG = G.copy()

        # find the ends of all the carbon tails
        tails = [tail[[len(list(nx.neighbors(G, n))) for n in tail].index(1)] for tail in ctails]
        logger.debug(f'tailn {tails}')
        if len(tails) == 2:  # this is a "standard" lipid with two fatty acid tails
            queryheads = [k for k, n in WG.nodes.items() if n['element'] != 'C' and n['element'] != 'O']
            logger.debug(f'queryheads {queryheads}')
            if len(queryheads) == 0:  # only non-carbons in head are O's
                queryheads = [k for k, n in WG.nodes.items() if n['element'] == 'O']
            maxpathlen = [-1] * len(cbonds)
            claimedhead = [''] * len(cbonds)
            for nc in queryheads:
                for i, b in enumerate(cbonds):
                    path = nx.shortest_path(WG, source=b[0], target=nc)
                    pathlen = len(path)
                    if pathlen > maxpathlen[i]:
                        maxpathlen[i] = pathlen
                        claimedhead[i] = nc
            logger.debug(f'claimedhead {claimedhead}')
            headcontenders = list(set(claimedhead))
            headvotes = np.array([claimedhead.count(x) for x in headcontenders])
            heads = [headcontenders[np.argmax(headvotes)]]
            logger.debug(f'heads {heads}')
        else:
            # this could be a wacky multiple-chain lipid with an extended head group
            # like cardiolipin or one of those bacterial glycolipids
            OG = WG.copy()
            for tail in ctails:
                for atom in tail:
                    OG.remove_node(atom)
            mins = len(OG) ** 2
            headatom = None
            for atom in OG:
                g = OG.copy()
                g.remove_node(atom)
                c = [len(x) for x in list(nx.connected_components(g))]
                while 1 in c:
                    c.remove(1)
                c = np.array(c)
                # logger.debug(f'testing atom {atom} -> {c}')
                if len(c) > 1:
                    s = c.std()
                    if s < mins:
                        # logger.debug(f'c {c} {atom} {s}<{mins}')
                        headatom = atom
                        mins = s
            assert headatom != None
            logger.debug(f'headatom {headatom}')
            heads = [headatom]

        self.annotation['heads'] = heads
        self.annotation['tails'] = tails
        self.annotation['shortest_paths'] = {head: {tail: len(nx.shortest_path(WG, source=head, target=tail)) for tail in tails} for head in heads}

    def to_file(self, f):
        """ 
        This method writes the residue record to a file object.
        It writes the title card with the residue name, charge, and synonym,
        followed by the atom, bond, internal coordinate, and delete atom records.
        If there are no atoms, bonds, ICs, or delete records, it writes an empty line for that section.
        
        Parameters
        ----------
        f : file-like object
            A file-like object to which the residue record will be written.
        """
        f.write(self.blockstring)

    def show_error(self, buf):
        """ 
        This method checks the error code of the residue record and writes an error message to the provided buffer if the error code is not 0.
        It provides specific error messages based on the error code value.
        
        Parameters
        ----------
        buf : function
            A function that takes a string as input and writes it to a buffer or log.
        """
        if self.error_code == -5:
            buf(f'{self.resname} references atoms in IC\'s that are not declared at ATOM\'s')

    def to_psfgen(self, W: PsfgenScripter, **kwargs):
        """ 
        Generate the psfgen script necessary to build the residue as as standalone molecule to generate its PSF/PDB.
        
        Parameters
        ----------
        W : PsfgenScripter
            A PsfgenScripter object to which the residue record will be written in PSFGEN format.
        **kwargs : dict
            Additional keyword arguments that can be passed to the method.
            
            - refic-idx: An index to select a specific internal coordinate (IC) to use as a reference for setting coordinates. Defaults to 0.
        
        Returns
        -------
        int
            Returns 0 if successful, or -1 if no valid reference atom could be found.
        """
        refic_idx = kwargs.get('refic-idx', 0)
        # find the first heavy atom that is atom-i in an IC w/o H that is not an improper
        refatom = None
        refic = None
        logger.debug(f'trying to find a reference atom from {len(self.atoms)} atoms in this resi')
        try:
            is_proper = [((not x.empty) and (x.dihedral_type == 'proper') and (not any([self.atomdict[a].element == 'H' for a in x.atoms]))) for x in self.IC]
        except:
            logger.warning(f'Could not parse IC\'s')
            return -1
        workingIC = list(compress(self.IC, is_proper))
        for ic in workingIC:
            logger.debug(f'{str(ic)}')
        if not workingIC:
            logger.warning(f'No valid IC for {self.resname}')
            return -1
        nWorkingIC = len(workingIC)
        logger.debug(f'{nWorkingIC}/{len(self.IC)} IC\'s with proper dihedrals and no hydrogens')
        refic = workingIC[refic_idx]
        refatom = self.atomdict[refic.atoms[1]]
        if refatom is None:
            logger.debug(f'Unable to identify a reference atom')
            return -1
        logger.debug(f'refIC is {str(refic)}')
        logger.debug(f"{[self.atomdict[a].element=='H' for a in refic.atoms]}")
        logger.debug(f'refatom is {refatom.name}')
        W.addline(r'segment A {')
        W.addline(r'    first none')
        W.addline(r'    last none')
        W.addline(f'    resid 1 {self.resname}')
        W.addline(r'}')
        # put reference atom at origin
        W.addline(f'psfset coord A 1 {refatom.name} '+r'{0 0 0}')
        b1a1, b1a2, b1l = refic.bonds[0].name1, refic.bonds[0].name2, refic.bonds[0].length
        partner = b1a1
        # partner is the first atom in the IC, put along x-axis
        W.addline(f'psfset coord A 1 {partner} '+r'{'+f'{b1l:.5f}'+r' 0 0}')
        a1a1, a1a2, a1a3, a1d = refic.angles[0].name1, refic.angles[0].name2, refic.angles[0].name3, refic.angles[0].degrees
        W.addline(f'set a [expr {a1d}*acos(-1.0)/180.0]')
        W.addline(f'set x [expr {b1l}*cos($a)]')
        W.addline(f'set y [expr {b1l}*sin($a)]')
        W.addline(f'psfset coord A 1 {a1a3} [list $x $y 0.0]')
        return 0

class CharmmResiList(UserList[CharmmResi]):
    """
    A list of CharmmResi objects.
    """
    def to_dict(self):
        return {r.resname: r for r in self}

    def from_dict(self, d: dict[str, CharmmResi]):
        self.clear()
        self.extend(d.values())

    @classmethod
    def from_blockstring_list(cls, blockstring_list: list[str], metadata: dict[str, str]):
        resi_list = cls()
        for block in blockstring_list:
            resi = CharmmResi.from_blockstring(block, metadata=metadata)
            resi_list.append(resi)
        return resi_list
    
    def to_resi_pres(self):
        resi = CharmmResiList(list(filter(lambda r: r.key=='RESI', self)))
        pres = CharmmResiList(list(filter(lambda r: r.key=='PRES', self)))
        return resi, pres

    def tally_masses(self, massdict: CharmmMassDict):
        for resi in self.data:
            resi.set_mass(massdict)

class CharmmResiDict(UserDict[CharmmResi]):
    """
    A dictionary mapping residue names to their corresponding CharmmResi objects.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @classmethod
    def from_blockstring_list(cls, blockstring_list: list[str], metadata: dict[str, str], resnames: list[str] = []):
        resi_dict = cls()
        for block in blockstring_list:
            resname = CharmmResi.name_from_blockstring(block)
            # logger.debug(f'Look-ahead found resname {resname}; include? {"yes" if len(resnames)==0 or resname in resnames else "no"}')
            if len(resnames) == 0 or resname in resnames:
                resi = CharmmResi.from_blockstring(block, metadata=metadata)
                resi_dict[resi.resname] = resi
        return resi_dict
    
    def to_resi_pres(self):
        resi = CharmmResiDict({k: v for k, v in self.items() if v.key == 'RESI'})
        pres = CharmmResiDict({k: v for k, v in self.items() if v.key == 'PRES'})
        return resi, pres

    def to_dict(self):
        return {k: v for k, v in self.items()}

    def from_dict(self, d: dict[str, CharmmResi]):
        self.clear()
        self.update(d)

    def tally_masses(self, massdict: CharmmMassDict):
        """
        Tally the masses of all residues in the dictionary.
        """
        for resi in self.values():
            resi.set_mass(massdict)

    def get_residue(self, name: str) -> CharmmResi | None:
        """
        Get a CharmmResi object by its name.

        Parameters
        ----------
        name : str
            The name of the residue to retrieve.

        Returns
        -------
        CharmmResi | None
            The corresponding CharmmResi object, or None if not found.
        """
        return self.get(name, None)