.. _example bpti2:

Example 2: BPTI Excluding the Phosphate Ion
-------------------------------------------

This is the same as Example 1, except we delete the phosphate ion.

.. literalinclude:: ../../../../pestifer/resources/examples/ex02/inputs/bpti2.yaml
    :language: yaml

Note the ``exclude`` subdirective under ``source``.  You remember how you can learn about it?  Using ``config-help``: 

.. code-block:: bash

  $ pestifer --no-banner config-help tasks psfgen source exclude

    exclude:
        Logical expressions involving atom attributess used in parallel to
        specify those atoms to be excluded

    All subattributes at the same level as 'exclude':

    base|tasks->psfgen->source
        biological_assembly
        transform_reserves
        remap_chainIDs
        reserialize
        model
        cif_residue_map_file
        include
        exclude
        sequence ->
        .. up
        ! quit
    pestifer-help:

The logical expressions are expected to by Pythonic, not VMD/Tcl-like.  In this case, we exclude any residue with the name "PO4".  You can use any atom attributes in the expression, such as ``resid``, ``chain``, ``name``, etc.  You can also combine multiple expressions using Python's logical operators ``and``, ``or``, and ``not``.  If you list multiple expressions under ``exclude``, they are combined with ``or`` logic.  So, for example, to exclude both phosphate ions and water molecules, you could use:

.. code-block:: yaml

    psfgen:
      source:
        exclude:
          - name == "PO4"
          - resname == "HOH"

Because the expressions can be compound, this would be equivalent to the single expression:

.. code-block:: yaml

    psfgen:
      source:
        exclude:
          - name == "PO4" or resname == "HOH"

We have also modified the ``solvate`` task to allow for a 0.154 M NaCl solution, which is a common salt concentration in biological systems.

.. _digression-validate-task:

Digression:  The Validate Task
==============================

This example also uses the `validate <_subs_runtasks_validate>`_ task, which is a useful way to check that your configuration file is doing what you expect by directly interrogating the PSF and PDB file of the current state.  This particular test validates the exclusion of the phosphate ion.  Other types of tests can check for presence or absence of other residues, interresidue bonds (disulfides and glycosylations), and more.

.. raw:: html

        <div class="autogen-footer">
            <p>Example author: Cameron F. Abrams&nbsp;&nbsp;&nbsp;Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
        </div>