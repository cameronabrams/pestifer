.. _subs_runtasks_psfgen_mods_grafts:

grafts
------
..
        # shortcode format: "target:source"
        # target format:    "C_RRR[-SSS]"
        #  - C chainID 
        #  - RRR resid+insertion of first residue in target alignment basis
        #  - SSS optional resid+insertion of last residue in target alignment basis
        #    (if not present, only RRR is used)
        # source format:     "pdbid,C_RRR[#SSS][-TTT]"
        #  - pdbid basename of pdb file or pdb id
        #  - C chainID in source pdb
        #  - RRR resid+insertion of first residue in source alignment basis
        #  - SSS optional resid+insertion of last residue in target alignment basis
        #    (if not present, only RRR is used)
        #  - TTT optional resid+insertion completing RRR-TTT range for entire graft source
        #    (if not present, only RRR is the entire graft source)

A ``grafts`` directive is used to specify grafting of a structure onto a the main structure.  The syntax for specifying a graft is as follows:

``<target>:<source>``

The ``<target>`` format is ``<chain>_<resi>[-<resi>]``, where ``<chain>`` is the chain ID, ``<resi>`` is the residue number (including an optionally appended insertion code) of the first residue in the target alignment basis, and ``[-<resi>]`` is an optional second residue number (including an optionally appended insertion code) of the last residue in the target alignment basis.  For example, ``A_12-15`` would graft the source structure on top of residue 12 of chain A by aligning congruent residues in the source onto residues 12 through 15.

The ``<source>`` format is ``<pdbid>,<chain>_<resi>[#<resi>][-<resi>]``, where ``<pdbid>`` is the basename of the PDB file or PDB ID, ``<chain>`` is the chain ID in the source PDB, ``<resi>`` is the residue number (including an optionally appended insertion code) of the first residue in the source alignment basis, ``[#<resi>]`` is an optional second residue number (including an optionally appended insertion code) of the last residue in the source alignment basis, and ``[-<resi>]`` is an optional third residue number (including an optionally appended insertion code) of the last residue in the graft source.  For example, ``1A2B,A_12#15-20`` would graft residues 12 through 20 from chain A of PDB 1A2B using residues 12 through 15 of chain A as the alignment basis onto the target.

A good example of using grafts to fully glycosylate a protein is :ref:`example 13`.


