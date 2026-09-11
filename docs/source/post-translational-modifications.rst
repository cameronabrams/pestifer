.. _post_translational_modifications:

Post-translational modifications
================================

Most post-translational modifications need no new parameters, no patches, and no configuration
beyond a single line.  CHARMM ships several hundred modified amino acids as **whole residues**
(``RESI`` entries), pestifer already classifies them as protein, and a ``mutations`` directive
swaps one in:

.. code-block:: yaml

    - psfgen:
        mods:
          mutations:
            - A:TYR,59,TYS      # tyrosine 59 -> sulfotyrosine

That is the whole mechanism.  ``TYS`` is not a patch and not a ligand; it is a tyrosine with a
sulfate, defined in the force field you already ship, and psfgen builds it and bonds it into the
chain like any other residue.

A modification that is *already in the input* needs even less.  A deposit carrying ``SEP`` --
declared by a ``MODRES`` record -- builds with no configuration at all, because the residue is
already classified as protein and nothing aliases it away.  Selenomethionine is the deliberate
exception: ``MSE`` **is** aliased to ``MET``, because it is a phasing substitution rather than
chemistry anyone wants to simulate.

Modified residues are reached by ``mutations``, not ``patches``
---------------------------------------------------------------

This is the most common mistake, and the two are not interchangeable:

- ``mutations`` replaces a whole residue with another whole residue.  Use it for anything CHARMM
  defines as a ``RESI`` -- which is nearly every modified amino acid.
- ``patches`` applies a CHARMM ``PRES``, which *edits* a residue in place.  Use it for the
  protonation-state patches (``ASPP``, ``GLUP``, ``LSN``...), for the phospho patches when you
  want control of charge state (below), and for linkages.

Handing a residue name to ``patches`` is now reported rather than silently ignored:

.. code-block:: text

    ERROR> unknown patch name(s) in this task's `patches`: 'SEP' is a CHARMM residue (RESI),
           not a patch (PRES); a whole modified residue is applied as a mutation, not a patch

What works today
----------------

Every residue below is a shipped ``RESI``, classified as protein, in a topology file pestifer
loads by default.  Each is one ``mutations`` line.

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Modification
     - Residues
   * - Sulfation
     - ``TYS`` (Tyr), ``OSE`` (Ser), ``ALS`` (Ser)
   * - Phosphorylation, His and Asp
     - ``NEP``, ``HIP`` (His), ``PHD`` (Asp)
   * - Lysine acylation
     - ``ALY`` acetyl, ``KCR`` crotonyl, ``PRK`` propanoyl, ``KHB`` hydroxybutyryl, ``KCX`` carboxy, ``MCL`` carboxyethyl
   * - Lysine and arginine methylation
     - ``MLZ`` mono, ``MLY`` di, ``M3L`` tri; ``AGM``, ``2MR``, ``DA2``, ``NMM`` (Arg)
   * - Histidine methylation
     - ``HIC``, ``MHS``, ``MHSP``
   * - Hydroxylation
     - ``HYP``, ``HZP`` (Pro), ``LYZ`` (Lys), ``3GL`` (Glu), ``AHB`` (Asn), ``BHD`` (Asp)
   * - Cysteine oxidation and adducts
     - ``CSO``, ``CSX``, ``CSU``, ``OCS``, ``CSS``, ``SNC`` nitroso, ``CME``, ``SMC``
   * - Other
     - ``CGU`` gamma-carboxyglutamate, ``CIR`` citrulline

This is a selection, not the whole set -- the force field carries several hundred.  To check any
residue, ask:

.. code-block:: bash

    $ pestifer show-resources resname TYS
    TYS
      topology: residue (RESI) defined in toppar_all36_prot_modify_res.str (segtype: protein)
                "O-sulfo-L-tyrosine"
      source:   native to the CHARMM release
      PDB repo: no coordinates in the built-in PDB repository

``segtype: protein`` and a topology file pestifer loads are what make a residue usable as a
mutation target.  Use the PDB chemical component ID -- CHARMM's names for modified amino acids
mostly match it.

Two things that will bite you
-----------------------------

**The minimize is not optional.**  A modification introduced by mutation has no coordinates in
the input.  psfgen creates its atoms from the topology and ``guesscoord`` places them, initially
at *placeholder* bond lengths of exactly 1.0 Å:

.. code-block:: text

    straight out of psfgen        after a minimize        CHARMM
      OG-P    1.00 Å                1.54 Å                 ~1.60
      P-O1P   1.00 Å                1.49 Å                 ~1.51

Put an ``md`` minimize immediately after the ``psfgen`` task, and put any ``validate`` test
*after* the minimize -- checking bond lengths, not just that the atom exists.

**Check the charge state.**  A residue name carries one protonation state, and it may not be the
one you want.  ``RESI SEP`` is the **monoanionic** phosphoserine, while the dianion dominates at
pH 7.  Where the charge state matters, use the patch route instead: build the ordinary residue
and apply ``SP1`` (mono) or ``SP2`` (di) for phosphoserine, ``THP1``/``THPB`` for phosphothreonine, and
``TP1``/``TP2`` for phosphotyrosine.

Phosphorylation
---------------

``SEP`` (phosphoserine) and ``TPO`` (phosphothreonine) are mutation targets like any other:

.. code-block:: yaml

    mutations:
      - A:SER,65,SEP      # the PINK1 site of ubiquitin

They live in ``toppar_all36_prot_na_combined.str``, which is not among the topology files
pestifer loads by default -- but it does not have to be listed.  The topology that defines a
mutation target is pulled in automatically, so nothing beyond the line above is required.
:ref:`Example 28 <example phosphoubiquitin>` is the worked case.

Phosphohistidine (``NEP``, ``HIP``) and aspartyl phosphate (``PHD``) are ordinary mutation
targets too; they are defined in a stream that is loaded anyway.

Phosphotyrosine (``PTR``) works too, and its parameters are worth a note.  ``RESI PTR`` types
the phenol oxygen as ``ON2B``, a nucleic-acid ester oxygen, and the angles it needs are the
phenol-phosphate set published in 1994 -- ``CA ON2b P``, ``ON3 P ON2b``, ``ON4 P ON2b``.  Those
are written with a lowercase ``b`` in the force field while the ``MASS`` record, and therefore
the PSF, says ``ON2B``.  CHARMM atom types are case-insensitive, so this is perfectly legal; it
was pestifer that matched them case-sensitively and silently dropped the records, which made
phosphotyrosine fail at the first dynamics step with ``UNABLE TO FIND ANGLE PARAMETERS``.  Fixed
-- but the shape is worth knowing if you ever meet a missing-parameter error for a residue whose
parameters you can see in the force field with your own eyes.

Limitations
-----------

**Lipidation and prenylation build, with one caveat.**  ``CYSP`` (S-palmitoyl), ``CYSF``
(farnesyl), ``CYSG`` (geranylgeranyl), ``CYSL``, ``GLYM`` (N-myristoyl-glycine) and ``LYSM``
(N-myristoyl-lysine) are mutation targets like any other:

.. code-block:: yaml

    mutations:
      - A:CYS,14,CYSP     # palmitoylate a cysteine

The acyl chain is built from the topology's internal coordinates and comes out in a proper
extended conformation -- a palmitoyl measures 19.1 Å from C1 to C16, against ~19 Å fully
extended, with no clash against the protein.

The caveat is **placement, in a membrane build**.  pestifer has no machinery that puts a
protein-attached acyl or prenyl tail *into the bilayer*; the chain is built relative to the
residue it hangs off, wherever that points.  For a soluble protein that is the right answer and
the tail is simply solvent-exposed.  For a membrane system, check where the tail ended up rather
than assuming it inserted -- a palmitoylated cytoplasmic tail whose chain is lying in water is
a system that will run and will be wrong.

**Glycosylation is a different mechanism.**  Sugars are separate residues joined to the protein
by a patch, not substituted for it; see :ref:`links <subs_buildtasks_psfgen_mods_links>` and
:ref:`grafts <subs_buildtasks_psfgen_mods_grafts>`, and examples 7, 8, 10, 15 and 17.

.. note::

   The segtype of a modified residue is derived from *which force-field file defines it*, not
   from its chemistry.  That is right for almost everything, and wrong for anything defined in a
   non-protein stream.  The six lipidated amino acids were derived as ``lipid`` for exactly this
   reason -- a palmitoylated cysteine could not be built into its own chain -- and are now
   curated as protein.  ``LYX``, a coenzyme-A lysine adduct, still derives as ``glycan``.  If a
   modified amino acid behaves as though it is not part of the chain, check its segtype first
   with ``pestifer show-resources resname``; the fix is a curated entry, which wins over the
   derivation.
