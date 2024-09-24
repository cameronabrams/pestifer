``exclude``
===========

Specifies any residues or atoms present in the PDB source to exclude from the system

Single-valued parameters:

  * ``chains``: Specify list of chain IDs to ignore; in PDB-format input, these are typically author-generated, while in mmCIF they are not.  You must use chain IDs that reference your chosen coordinate file input format.

  * ``resnames``: Specify list of resnames to ignore; good for excluding waters or ions or small molecules.  Resnames must reference your chosen coordinate input file.

  * ``resseqnums``: Specify list of resseqnums to ignore



