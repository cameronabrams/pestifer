# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
On-the-fly CGenFF parameterization for ligands not covered by the
standard CHARMM topologies.

The pipeline:
  pestifer Residue (HETATM)
    -> protonate at target pH (RDKit + Dimorphite-DL)
    -> write mol2 (Open Babel, format conversion only)
    -> shell out to the ``cgenff`` binary
    -> parse the .str output
    -> remap CGenFF-invented heavy-atom names back to the PDB names
    -> inject into the psfgen session as an extra topology stream

This subpackage builds that pipeline piece by piece.
"""

from .protonation import LigandProtonationError, protonate_ligand
from .rcsb import RCSBLookupError, fetch_ligand_smiles

__all__ = [
    "LigandProtonationError",
    "protonate_ligand",
    "RCSBLookupError",
    "fetch_ligand_smiles",
]
