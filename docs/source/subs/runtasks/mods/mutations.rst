.. _subs_runtasks_psfgen_mods_mutations:

mutations
---------
..
        ### shortcode format: c:nnn,rrr,mmm or c:A###B
        ### c: chainID
        ### nnn: old resname, 1-byte rescode allowed
        ### rrr: resseqnum
        ### mmm: new resname, 1-byte rescode allowed

A ``mutations`` directive is used to specify point mutations, which are changes in a single amino acid residue.  The syntax for specifying a mutation is as follows:

``<chain>:<oldresi>,<seqnum>,<newresi>``

``<chain>`` is the chain ID, ``<oldresi>`` is the old residue name, ``<seqnum>`` is the residue number, and ``<newresi>`` is the new residue name.  You can use either 1-byte or 3-byte residue names.  For example, both ``A:ALA,12,GLU`` and ``A:A,12,E`` would change the alanine at position 12 in chain A to glutamic acid.