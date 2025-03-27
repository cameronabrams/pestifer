.. _subs_runtasks_psfgen_mods_ssbonds:

ssbonds and ssbondsdelete
-------------------------

ssbonds
^^^^^^^
..
        # shortcode format: C_RRR-D_SSS
        # C, D chainIDs
        # RRR, SSS resids

An ``ssbonds`` directive is used to specify disulfide bonds between two residues.  The syntax for specifying a disulfide bond is as follows:

``<chain1>_<resi1>-<chain2>_<resi2>``

``<chain1>`` and ``<chain2>`` are the chain IDs, and ``<resi1>`` and ``<resi2>`` are the residue numbers of the two residues that form the disulfide bond.  For example, ``A_12-B_15`` would create a disulfide bond between residue 12 in chain A and residue 15 in chain B.  You can use either 1-byte or 3-byte residue names.

ssbondsdelete
^^^^^^^^^^^^^

An ``ssbondsdelete`` directive is used to specify the deletion of a disulfide bond between two residues.  The syntax for specifying a disulfide bond deletion is as follows:

``<chain1>_<resi1>-<chain2>_<resi2>``

``<chain1>`` and ``<chain2>`` are the chain IDs, and ``<resi1>`` and ``<resi2>`` are the residue numbers of the two residues that form the disulfide bond.  For example, ``A_12-B_15`` would delete the disulfide bond between residue 12 in chain A and residue 15 in chain B.  You can use either 1-byte or 3-byte residue names.
