# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definea the :class:`PSFAtom` and :class:`PSFAtomList` classes for handling atoms in PSF topology files.
These classes are used to represent individual atoms and collections of atoms in a PSF file.
The :class:`PSFAtom` class represents a single atom, while the :class:`PSFAtomList` class is a collection of such atoms.
These classes provide methods for parsing atom lines and managing atom properties."""

import logging
import networkx as nx

from .psftopoelement import PSFTopoElement,PSFTopoElementList

from ..core.labels import Labels
from ..core.stringthings import split_ri

logger=logging.getLogger(__name__)

class PSFAtom(PSFTopoElement):
    """
    A class representing an atom in a PSF topology file.
    An atom is defined by its serial number, chain ID, residue sequence number, insertion code, residue name, atom name, type, charge, atomic weight, and segment type.
    The class provides methods for parsing atom lines, checking if the atom is a hydrogen atom,
    adding ligands, and checking if the atom is part of a peptide bond.

    Parameters
    ----------
    atomline : str
        A string representing a line from a PSF file that contains information about the atom.
        The line should contain the following fields in order: serial number, chain ID, residue sequence number, insertion code, residue name, atom name, type, charge, and atomic weight.
        The line is expected to have exactly 9 tokens when split by whitespace.
    """
    def __init__(self,atomline):
        tokens=[x.strip() for x in atomline.split()]
        assert len(tokens)==9,f'cannot parse psf atomline: {atomline}'
        idx=int(tokens[0])
        super().__init__([idx])
        r,i=split_ri(tokens[2])
        self.serial=idx
        self.chainID=tokens[1]
        self.resseqnum=r
        self.insertion=i
        self.resname=tokens[3]
        self.name=tokens[4]
        self.type=tokens[5]
        self.charge=float(tokens[6])
        self.atomicwt=float(tokens[7])
        self.segtype=Labels.segtype_of_resname[self.resname]
        self.ligands=[]
    
    def __hash__(self):
        return self.serial
    
    def __eq__(self,other):
        return self.serial==other.serial

    def __str__(self):
        return f'{self.chainID}_{self.resname}_{self.resseqnum}{self.insertion}-{self.name}({self.serial})'

    def isH(self):
        """
        Check if the atom is a hydrogen atom.
        Returns True if the atom's type is 'H', otherwise False.
        """
        return self.atomicwt<1.1
    
    def add_ligand(self,other):
        """
        Add a ligand to the atom's list of ligands.

        Parameters
        ----------
        other : PSFAtom
            The ligand atom to add.
        """
        if not other in self.ligands:
            self.ligands.append(other)

    def is_pep(self,other):
        """
        Check if the atom is part of a peptide bond with another atom.
        A peptide bond is defined as a bond between a carbon atom (C) and a nitrogen atom (NH1).

        Parameters
        ----------
        other : PSFAtom
            The other atom to check for a peptide bond.
        """
        return (self.type=='C' and other.type=='NH1') or (self.type=='NH1' and other.type=='C')

class PSFAtomList(PSFTopoElementList):
    """
    A class representing a list of :class:`PSFAtom` objects.
    This class inherits from :class:`PSFTopoElementList <.psftopoelement.PSFTopoElementList>` and provides methods for managing a collection of atoms in a PSF topology file.

    Parameters
    ----------
    atoms : list of PSFAtom
        A list of `PSFAtom` objects representing the atoms in the PSF file.

    """
    def __init__(self,atoms):
        super().__init__(atoms)

    def graph(self):
        """
        Create a :class:`networkx.Graph` representation of the atom list, where atoms are nodes and bonds are edges.
        This method constructs a networkx graph where each atom is a node and edges are created between atoms and their ligands.
        
        Returns
        -------
        networkx.Graph
            A networkx graph representing the atom list, with atoms as nodes and bonds as edges
        """
        g=nx.Graph()
        for a in self:
            for l in a.ligands:
                g.add_edge(a,l)
        logger.debug(f'atomlist graph has {len(g)} nodes')
        return g
