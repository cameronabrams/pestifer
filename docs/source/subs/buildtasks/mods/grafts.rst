.. _subs_buildtasks_psfgen_mods_grafts:

grafts
------

A ``grafts`` directive is used to specify grafting of a structure onto the main structure.  Each graft is expressed as a shortcode string:

``<target>:<source>``

Target specification
++++++++++++++++++++

The ``<target>`` has the format::

    <chain>_<R1>[#<R2>[#<R3>...]]

where ``<chain>`` is the chain ID in the base structure and ``<R1>``, ``<R2>``, etc. are residue numbers (with optional insertion codes) that serve as the **alignment basis** — the residues in the base structure onto which the congruent source residues will be superimposed.  Multiple residues are separated by ``#``.

Examples:

- ``A_512`` — single-residue target: align to residue 512 of chain A
- ``G_542#543`` — two-residue target: align to residues 542 and 543 of chain G
- ``G_572#573#574`` — three-residue target: align to residues 572, 573, and 574 of chain G

Source specification
++++++++++++++++++++

The ``<source>`` has the format::

    <pdbid>,<chain>_<S1>[#<S2>[#<S3>...]][-<E>]

where ``<pdbid>`` is the basename of the PDB file or PDB ID, ``<chain>`` is the chain ID in the source PDB, ``<S1>``, ``<S2>``, etc. are the **source index residues** used for alignment (one per target residue, in the same order), and the optional ``[-<E>]`` limits the graft to source residues up to and including ``<E>``.

Rules:

- The number of source index residues must equal the number of target residues (N), **or** N+1, in which case the last one is interpreted as ``source_end``.
- Source residues that are not index residues (and are within the ``source_end`` limit if given) are the **donor residues** that become new atoms in the built structure.

Examples:

- ``4byh,C_1-10`` — align to source residue 1, donate residues 2–10
- ``4b7i,C_1#2-8`` — align to source residues 1 and 2, donate residues 3–8
- ``4b7i,C_1#2#3-8`` — align to source residues 1, 2, and 3, donate residues 4–8

N-point alignment and dihedral pre-conditioning
++++++++++++++++++++++++++++++++++++++++++++++++

When more than one alignment residue is specified, pestifer performs **N-point alignment** with **dihedral pre-conditioning**.

This is particularly useful for glycan grafting onto structures that already carry a partial glycan stub at the target site.  For example, if the base structure already has a two-sugar NAG–NAG stub at an N-glycosylation site, a two-residue target (``G_542#543``) aligned to congruent two-residue source resids (``4b7i,C_1#2-8``) will produce a much better placement of the donated sugars than a one-residue alignment.

**Pre-conditioning** runs before the alignment transformation is computed.  For each glycosidic bond between consecutive source index residues (processed innermost-to-outermost):

1. The target's phi and psi angles (and omega for 1→6 and 2→6 linkages) are measured in the base structure.
2. The corresponding angles are measured in the source PDB.
3. The difference (target − source) is used to rotate the source glycan around the appropriate bond axis, moving the outer (donor) side while keeping the innermost index residue fixed.

After pre-conditioning, the standard ``measure fit`` least-squares alignment is applied to all N index residues simultaneously, and the donor residues are transformed accordingly.

Full example
++++++++++++

A two-residue graft that extends the NAG–NAG stub at residue 542–543 of chain G with residues 3–8 from PDB 4b7i chain C::

    - G_542#543:4b7i,C_1#2-8

And a three-residue graft for a site where the base already has a three-sugar stub (NAG–NAG–BMA at 572–573–574)::

    - G_572#573#574:4b7i,C_1#2#3-8

A good example of using grafts to fully glycosylate a protein is :ref:`example sars cov2 spike ba2`.
