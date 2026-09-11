.. _subs_buildtasks_psfgen_mods_mutations:

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

Mutating to a modified residue
++++++++++++++++++++++++++++++

``<newresi>`` does not have to be one of the twenty standard amino acids.  CHARMM ships several
hundred **modified** amino acids as whole residues -- phosphoserine ``SEP``, sulfotyrosine
``TYS``, acetyl-lysine ``ALY``, methyl-lysine ``MLZ``, hydroxyproline ``HYP`` and so on -- and a
mutation is how one is installed:

.. code-block:: yaml

    mutations:
      - A:TYR,59,TYS      # tyrosine 59 -> sulfotyrosine

This is the route for post-translational modifications; ``patches`` is not (it applies a CHARMM
``PRES``, which edits a residue rather than replacing it).  The modification has no coordinates
in the input, so a ``md`` minimize must follow the ``psfgen`` task.  See
:ref:`Post-translational modifications <post_translational_modifications>` for what is available,
the charge-state caveat, and the one family that needs an extra topology stream.

