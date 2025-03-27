.. _subs_runtasks_psfgen_mods_links:

links
-----

A ``links`` directive is used to specify the addition of a covalent bond between two atoms.  The syntax for specifying a link is as follows:

``<chain1>_<resi1>_<atom1>-<chain2>_<resi2>_<atom2>``

``<chain1>`` and ``<chain2>`` are the chain IDs, ``<resi1>`` and ``<resi2>`` are the residue numbers of the two residues that form the link, and ``<atom1>`` and ``<atom2>`` are the atom names of the two atoms that form the link.  For example, ``A-12_O-C_15_N`` would create a covalent bond between the oxygen atom of residue 12 in chain A and the nitrogen atom of residue 15 in chain C.  You can use either 1-byte or 3-byte residue names.

Links are not normally needed in a ``mods`` directive of a ``psfgen`` task.  However, if you are using a non-published pdb file that does not include metadata indicating non-protein bonds (e.g., a covalent bond between a glycan monomer and a protein), you can use the ``links`` directive to specify each required bond.  The ``links`` directive is not needed for standard protein structures.
