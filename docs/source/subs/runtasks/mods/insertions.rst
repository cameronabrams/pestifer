.. _subs_runtasks_psfgen_mods_insertions:

insertions
----------
..
        # C,R,SSS
        # C: chain ID
        # R: resid+insertioncode
        # SSS: one-letter amino acid sequence to insert

An ``insertions`` directive is used to specify insertions of a single residue.  The syntax for specifying an insertion is as follows:

``<chain>,<resid>,<resname>``

``<chain>`` is the chain ID, ``<resid>`` is the residue number (including an optionally appended insertion code) immediately after which the inserted residue is to be placed, and ``<resname>`` is the one-letter amino acid sequence to insert.  The inserted residue acquires the same integer residue number as the residue indicated and an incremented insertion code (usually ``A`` if the reference residue itself does not already have an insertion code). For example, ``A,12,A`` would insert an alanine at position 12A in chain A .  You can use either 1-byte or 3-byte residue names.  For example, both ``A,12,A`` and ``A,12,ALA`` would insert an alanine at position 12A in chain A.

