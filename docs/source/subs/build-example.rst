.. _sub_build_example:

build-example
-------------

There are 26 example systems; to run number four, for example:

.. code-block:: bash

   $ pestifer build-example 4

.. note::

   ``run-example`` is accepted as a synonym for ``build-example`` for backwards compatibility.

(Best to do that in a clean directory.)  

You can use ``pestifer show-resources`` to list the examples:

.. code-block:: bash

    $ pestifer show-resources examples

      Examples:

           ID      DBID  Name                            Title
            1      6pti  bpti1                           Bovine Pancreatic Trypsin Inhibitor (BPTI)
            2      6pti  bpti2                           BPTI with phosphate ion excluded in salty solution
            3      6pti  bpti3                           BPTI, no phosphate, some random mutations plus deletion of one disulfide
            4      6pti  bpti4                           BPTI, no phosphate, introducing a disulfide via mutations
            5      1f7a  hiv-protease                    HIV-1 protease dimer
            6      1fas  green-mamba-toxin               Fasciculin 1, an Anti-Acetylcholinesterase Toxin from Green Mamba Snake Venom
            7      4zmj  hiv-sosip-env-ectodomain1       Closed, Unliganded HIV-1 BG505 Env SOSIP-664 Trimer
            8      4tvp  hiv-sosip-env-ectodomain2       Closed, PGT122/35O22-Liganded HIV-1 BG505 Env SOSIP.664 Trimer (Fabs removed)
            9      8fad  hiv-ad8-env-ectodomain          HIV-1 Env Trimer 8fad
           10      8fae  hiv-ae2-env-ectodomain          HIV-1 Env Trimer 8fae
           11      7txd  hiv-sosip-env-ectodomain3       HIV-1 Env Trimer 7txd, substitutions, no ligands
           12      5vn3  hiv-sosip-env-ectodomain4       HIV-1 Env Trimer 5vn3, excluding sCD4 and Fab chains, with Gly3 stubs replacing missing v1/2
           13      2ins  insulin-hexamer                 hexameric insulin
           14      4zxb  insulin-receptor-ectodomain     insulin receptor ectodomain, Fabs removed
           15      7xix  sars-cov2-S-BA2                 BA.2 SARS-CoV-2 Spike 7xix, fully glycosylated using grafts, and cleaved
           16      6e8w  hiv-mpertm3-membrane1           HIV-1 gp41 MPER-TM trimer 6e8w embedded in DMPC bilayer (grid packer)
           17      6e8w  hiv-mpertm3-membrane2           HIV-1 gp41 MPER-TM trimer 6e8w embedded in model viral bilayer (grid packer)
           18      5fkw  ecoli-polymerase                E. Coli replicative DNA polymerase complex bound to a primer-template DNA
           19      1mob  sperm-whale-myoglobin           Sperm whale myoglobin
           20    P22033  methylmalonyl-coa-mutase        Mitochondrial methylmalonyl-CoA mutase (Alphafold P22033)
           21      1aon  groel-groes-adp                 Asymmetric GroEL/GroES Chaperonin Complex
           22      2bgj  ferredoxin-fad                  Ferredoxin-NADP(H) Reductase from Rhodobacter capsulatus
           23      1scd  subtilisin-dmso                 Subtilisin Carlsberg in DMSO
           24      1scd  subtilisin-acetone              Subtilisin Carlsberg in acetone
           25       N/A  hiv-env-cd4-17b-liganded        Open HIV-1 B41 SOSIP Env trimer with sCD4 and completed 17b Fabs (5vn3)
           26      1scd  subtilisin-acetonitrile         Subtilisin Carlsberg in acetonitrile

.. _sub_fetch_example:

fetch-example
-------------

This subcommand is like ``build-example``, except it only copies the YAML input file needed to run the example to the CWD.  It can be edited, if desired, and then run using the ``build`` subcommand.

