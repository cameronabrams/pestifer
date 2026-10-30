.. _subs_buildtasks_psfgen_mods_substitutions:

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

``<chain>`` is the chain ID, ``<startresi>`` and ``<endresi>`` are the first and last residue numbers, respectively, of the subchain to be substitued, and ``<sequence>`` is the one-letter amino acid sequence to replace this subchain with.  For example, ``A:12-15,ACDE`` would replace residues 12 through 15 in chain A with alanine, cysteine, aspartic acid, and glutamic acid.  Letters go on in sequence order, whether a residue in the range is resolved (it is mutated in place, keeping its backbone) or missing (it is built as the new residue).  A sequence shorter than the range deletes the residues left over; one longer than the range does not add residues -- the extra letters are reported and dropped -- so use ``insertions`` for that.

