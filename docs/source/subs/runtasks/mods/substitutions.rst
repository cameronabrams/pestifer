.. _subs_runtasks_psfgen_mods_subtitutions:

substitutions
-------------
..
        # shortcode format: C:nnn-ccc,abcdef
        # C -- chainID
        # nnn -- N-terminal resid of sequence to be replaced
        # ccc -- C-terminal resid of sequence to be replaced
        # abcdef -- the 1-byte rescode sequence to be substituted

A ``substitutions`` directive is used to specify substitutions of a sequence of residues.  The syntax for specifying a substitution is as follows:

``<chain>:<startresi>-<endresi>,<sequence>``

``<chain>`` is the chain ID, ``<startresi>`` and ``<endresi>`` are the first and last residue numbers, respectively, of the subchain to be substitued, and ``<sequence>`` is the one-letter amino acid sequence to replace this subchain with.  For example, ``A:12-15,ACDEFG`` would replace residues 12 through 15 in chain A with alanine, cysteine, aspartic acid, glutamic acid, phenylalanine, and glycine.

