Example 2: BPTI with Phosphate Ion Deleted
------------------------------------------
This is the same as Example 1, except we delete the phosphate ion.

.. code-block:: yaml

  title: BPTI with phosphate ion excluded
  tasks:
    - psfgen:
        source:
          id: 6pti
          exclude:
            resnames:
              - PO4
    - md:
        ensemble: minimize
    - solvate:
    - md:
        ensemble: minimize
    - md:
        ensemble: NVT
    - md:
        ensemble: NPT
        nsteps: 200
    - md:
        ensemble: NPT
        nsteps: 400
    - md:
        ensemble: NPT
        nsteps: 800
    - md:
        ensemble: NPT
        nsteps: 1600
    - terminate:
        basename: my_6pti
        package:
          ensemble: NPT
          basename: prod_6pti

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
