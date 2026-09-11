.. _example phosphoubiquitin:

Example 28: Ser65-phosphorylated ubiquitin
------------------------------------------

`PDB ID 1ubq <https://www.rcsb.org/structure/1UBQ>`_ is ordinary human ubiquitin.  Nothing in it is phosphorylated.  This example installs a post-translational modification that is *not* in the input: Ser65 is converted to phosphoserine by mutating it to ``SEP``, the CHARMM residue for phosphoserine.

Ser65 is the PINK1 site.  Its phosphorylation activates parkin and is the committing step of mitophagy, so it is a modification one frequently wants to simulate on a structure that does not carry it.

.. code-block:: yaml

    mods:
      mutations:
        - A:SER,65,SEP

A modified residue is reached by ``mutations``, not by ``patches``: CHARMM ships several hundred modified amino acids -- ``SEP``, ``TPO``, ``PTR``, ``TYS``, ``MLZ`` and the rest -- as whole ``RESI`` residues rather than as patches, and a mutation is how a whole residue is swapped in.  (``patches`` is for true ``PRES`` entries; giving it a residue name is now reported rather than silently ignored.)

Two things about this build are worth reading before adapting it.

**The topology that defines the residue is pulled in automatically.**  ``SEP`` lives in ``toppar_all36_prot_na_combined.str``, which is not among the files pestifer loads by default; the topology defining a mutation target is added for you, so this config lists no force-field files at all.  (``TPO`` shares that file and works the same way.  ``PTR`` does not work in this CHARMM release -- see the warning in :ref:`Post-translational modifications <post_translational_modifications>`.)

**The minimize is not optional.**  The phosphate has no coordinates in the input; psfgen creates those atoms from the topology and ``guesscoord`` places them, initially at placeholder bond lengths of exactly 1.0 Å:

.. code-block:: text

    straight out of psfgen        after the minimize
      OG-P    1.00 A                OG-P    1.54 A     (CHARMM ~1.60)
      P-O1P   1.00 A                P-O1P   1.49 A     (CHARMM ~1.51)
      P-O2P   1.00 A                P-O2P   1.46 A     (CHARMM ~1.51)
      P-O3P   1.00 A                P-O3P   1.59 A     (CHARMM ~1.58)

A build that skipped straight from the mutation to dynamics would start from a collapsed phosphate.  The ``validate`` task therefore runs *after* the minimize, and checks that the phosphorus exists rather than assuming the mutation implies it.

.. note::

   ``RESI SEP`` is the **monoanionic** phosphoserine.  The dianion dominates at pH 7, so a study in which the charge state matters should instead build an ordinary ``SER`` and apply the ``SP2`` patch (``SP1`` gives the monoanion explicitly).  The same choice exists for phosphothreonine (``THP1``/``THPB``) and phosphotyrosine (``TP1``/``TP2``).  The mutation route used here matches the residue a depositor would have annotated; the patch route gives control of protonation.

A modification that is *already* in the input needs none of this.  A deposit carrying ``SEP`` -- declared by a ``MODRES`` record -- builds with no configuration at all, because the residue is already classified as protein and nothing aliases it away.  Selenomethionine is the deliberate exception: ``MSE`` **is** aliased to ``MET``, since it is a phasing substitution rather than chemistry to preserve.

Most other modifications need no extra stream either: sulfotyrosine, the lysine acyl and methyl
ladders, hydroxyproline and the cysteine oxidation series are all reachable with a single
``mutations`` line and nothing else.  :ref:`Post-translational modifications
<post_translational_modifications>` lists what is available and why phosphoserine is one of the
few that needs the configuration above.

.. literalinclude:: ../../../../pestifer/resources/examples/28/inputs/phosphoubiquitin.yaml
    :language: yaml

.. task-table:: ../../../../pestifer/resources/examples/28/inputs/phosphoubiquitin.yaml

.. raw:: html

    <div class="autogen-footer">
        <p>Example author: Cameron F. Abrams &nbsp;&nbsp;&nbsp; Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
    </div>
