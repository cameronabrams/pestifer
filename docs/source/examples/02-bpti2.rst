.. _example 2:

Example 2: BPTI with Phosphate Ion Deleted
------------------------------------------

This is the same as Example 1, except we delete the phosphate ion.

.. literalinclude:: ../../../pestifer/resources/examples/02-bpti-exclude-phosphate.yaml
    :language: yaml

Note the ``exclude`` subdirective under ``source``.  You remember how you can learn about it?  Using ``config-help``: 

.. code-block:: console

  $ pestifer --no-banner config-help tasks psfgen source exclude
  Help on user-provided configuration file format

  exclude:
      Specifies any residues or atoms present in the PDB source to exclude
        from the system

  base|tasks->psfgen->source->exclude
      chains
      resnames
      resseqnums
      .. up
      ! quit
  pestifer-help: resnames

  resnames:
      Specify list of resnames to ignore; good for excluding waters or ions
        or small molecules.  Resnames must reference your chosen
        coordinate input file.

  All subdirectives at the same level as 'resnames':

  base|tasks->psfgen->source->exclude
      chains
      resnames
      resseqnums
      .. up
      ! quit
  pestifer-help:

Each of ``chains`` and ``resnames`` are lists, and in the configuration file above, we have a single-element list for ``resnames`` that indicates the resname ``PO4``, which is how the phosphate ion is labelled in the original PDB file.

We have also modified the ``solvate`` task to allow for a 0.154 M NaCl solution, which is a common salt concentration in biological systems.
