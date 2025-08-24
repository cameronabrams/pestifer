# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Defines the :class:`PSFAtom` and :class:`PSFAtomList` classes for handling atoms in PSF topology files.
These classes are used to represent individual atoms and collections of atoms in a PSF file.
The :class:`PSFAtom` class represents a single atom, while the :class:`PSFAtomList` class is a collection of such atoms.
These classes provide methods for parsing atom lines and managing atom properties.
"""

from __future__ import annotations
import logging
import networkx as nx

from ..core.baseobj import BaseObj, BaseObjList
from pydantic import Field

from ..core.labels import Labels
from ..objs.resid import ResID

logger = logging.getLogger(__name__)

class PSFAtom(BaseObj):
    """
    A class representing an atom in a PSF topology file.
    An atom is defined by its serial number, chain ID, residue sequence number, insertion code, residue name, atom name, type, charge, atomic weight, and segment type.
    The class provides methods for parsing atom lines, checking if the atom is a hydrogen atom,
    adding ligands, and checking if the atom is part of a peptide bond.
    """

    _required_fields = { 'serial', 'resid', 'segname', 'resname', 'atomname', 'atomtype', 'charge', 'atomicwt', 'segtype' }
    _optional_fields = { 'ligands' }

    serial: int = Field(..., description="The serial number of the atom.")
    resid: ResID = Field(..., description="The residue ID of the atom, including residue number and insertion code.")
    segname: str = Field(..., description="The segment name of the atom.")
    resname: str = Field(..., description="The residue name of the atom.")
    atomname: str = Field(..., description="The name of the atom.")
    atomtype: str = Field(..., description="The type of the atom.")
    charge: float = Field(..., description="The charge of the atom.")
    atomicwt: float = Field(..., description="The atomic weight of the atom.")
    segtype: str = Field(..., description="The segment type of the atom.")
    ligands: list[PSFAtom] = Field(default_factory=list, description="A list of ligands associated with the atom.")

    @classmethod
    def _adapt(cls, *args, **kwargs) -> dict:
        if args and isinstance(args[0], str):
            atomline = args[0]
            return cls._from_psf_atomline(atomline)
        return super()._adapt(*args, **kwargs)

    @staticmethod
    def _from_psf_atomline(atomline: str):
        tokens = [x.strip() for x in atomline.split()]
        assert len(tokens) == 9, f'Cannot parse psf atomline: {atomline}'
        return PSFAtom(
            serial=int(tokens[0]),
            resid=ResID(tokens[2]),
            segname=tokens[1],
            resname=tokens[3],
            atomname=tokens[4],
            atomtype=tokens[5],
            charge=float(tokens[6]),
            atomicwt=float(tokens[7]),
            segtype=Labels.segtype_of_resname[tokens[3]],
        )
    
    def __hash__(self):
        return self.serial

    def __eq__(self, other: PSFAtom):
        return self.serial == other.serial

    def __str__(self):
        return f'{self.segname}_{self.resname}_{self.resid.resid}-{self.atomname}({self.serial})'

    def isH(self):
        """
        Check if the atom is a hydrogen atom.
        Returns True if the atom's type is 'H', otherwise False.
        """
        return self.atomicwt < 1.1

    def add_ligand(self, other: PSFAtom):
        """
        Add a ligand to the atom's list of ligands.

        Parameters
        ----------
        other : PSFAtom
            The ligand atom to add.
        """
        if not other in self.ligands:
            self.ligands.append(other)

    def is_pep(self, other: PSFAtom):
        """
        Check if the atom is part of a peptide bond with another atom.
        A peptide bond is defined as a bond between a carbon atom (C) and a nitrogen atom (NH1).

        Parameters
        ----------
        other : PSFAtom
            The other atom to check for a peptide bond.
        """
        return (self.atomtype == 'C' and other.atomtype == 'NH1') or (self.atomtype == 'NH1' and other.atomtype == 'C')

class PSFAtomList(BaseObjList[PSFAtom]):
    """
    A class representing a list of :class:`PSFAtom` objects.
    This class inherits from :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` and provides methods for managing a collection of atoms in a PSF topology file.
    """

    def describe(self) -> str:
        return f'PSFAtomList with {len(self.data)} atoms'
    
    def graph(self):
        """
        Create a :class:`networkx.Graph` representation of the atom list, where atoms are nodes and bonds are edges.
        This method constructs a networkx graph where each atom is a node and edges are created between atoms and their ligands.
        
        Returns
        -------
        networkx.Graph
            A networkx graph representing the atom list, with atoms as nodes and bonds as edges
        """
        g = nx.Graph()
        for a in self.data:
            for l in a.ligands:
                g.add_edge(a, l)
        logger.debug(f'atomlist graph has {len(g)} nodes')
        return g
