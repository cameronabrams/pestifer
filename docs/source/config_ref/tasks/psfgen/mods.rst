.. _config_ref tasks psfgen mods:

``mods``
========

Specifies modifications to the sequence or structure in the PDB file

Single-valued parameters:

  * ``mutations``: Single-residue mutations, expressed as shortcodes; e.g., C:F32A for mutation of Phe 32 on chain C to an Alanine

  * ``ssbonds``: additional ssbonds, expressed as shortcodes; e.g., C_154-A_298 to add an ssbond between resids 154 of chain C and 298 of chain A.  Note that the indicated residues must by CYS's.  If they are not, you should include each in a mutation shortcode in the mutations list

  * ``ssbondsdelete``: existing ssbonds to break, expressed as shortcodes; e.g., C_154-A_298 will break the existing ssbond between resids 154 of chain C and 298 of chain A so that each is a reduced CYS.

  * ``links``: additional generic bonds, expressed as shortcodes; e.g., C_154_N-A_298_M to add a link between atom N of resid 154 of chain C and atom M of resid 298 of chain A.

  * ``deletions``: residue ranges to delete, expressed as shortcodes; e.g., C:99-101 deletes residues 99 to 101 (inclusive) of chain C

  * ``substitutions``: residue ranges to replace with a new sequence; e.g., C:99-101,GG replaces residues 99, 100, and 101 of chain C with two glycines

  * ``insertions``: sequence of residues to insert C-terminal to a specified residue in a specified chain; e.g., C,123,GRETA inserts the residues Gly123A, Arg123B, Glu123C, Thr123D, and Ala123E immediately C-terminal to residue 123 of chain C; if the resid includes a terminal "+", then inserted resids are incremented as integers

  * ``crotations``: dihedral angle rotations to reposition modelled-in protein chains or glycans

  * ``Cfusions``: fuse residues from named residue range of named protein chain of named input PDB coordinate file to C-terminus of named segment of base molecule

  * ``grafts``: graft residues from named residue range of named chain of named input PDB coordinate file onto target residue of base molecule



