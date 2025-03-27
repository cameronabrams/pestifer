.. _subs_runtasks_psfgen_mods_deletions:

deletions
---------
..
        # shortcode format: C:nnn-ccc
        # C -- chainID
        # nnn -- N-terminal resid of sequence to be deleted
        # ccc -- C-terminal resid of sequence to be deleted

A ``deletions`` directive is used to specify deletions of a sequence of residues.  The syntax for specifying a deletion is as follows:

``<chain>:<startresi>-<endresi>``

``<chain>`` is the chain ID, ``<startresi>`` is the first residue number in the deletion, and ``<endresi>`` is the last residue number in the deletion.  For example, ``A:12-15`` would delete residues 12 through 15 in chain A.  You can use either 1-byte or 3-byte residue names.  For example, both ``A:12-15`` and ``A:A12-15`` would delete residues 12 through 15 in chain A.

