.. _sub_run_example:

run-example
-----------

There are 22 example systems; to run number four, for example:

.. code-block:: bash
   
   $ pestifer run-example 4

(Best to do that in a clean directory.)  

You can use ``pestifer show-resources`` to list the examples:

.. code-block:: bash

    $ pestifer show-resources examples

      Examples:

      Index        ID  Name                            Description
          1      6pti  bpti1                           Bovine Pancreatic Trypsin Inhibitor (BPTI)
          2      6pti  bpti2                           BPTI Excluding the Phosphate Ion
          3      6pti  bpti3                           BPTI with a Mutated-out Disulfide Bond
          4      6pti  bpti4                           BPTI with a Mutated-in Disulfide Bond
          5      1f7a  hiv-protease                    HIV Protease with Patches to Protonate Aspartates
          6      1fas  green-mamba-toxin               Green Mamba Toxin at pH 7.0
          7      4zmj  hiv-sosip-env-ectodomain1       Closed, Unliganded HIV-1 BG505 Env SOSIP.664 Trimer
          8      4tvp  hiv-sosip-env-ectodomain2       Closed, PGT122/35O22-Liganded HIV-1 BG505 Env SOSIP.664 Trimer (ligands removed)
          9      8fad  hiv-ad8-env-ectodomain          Cleaved, Asymmetric HIV-1 AD8 Env Ectodomain Trimer
         10      8fae  hiv-ae2-env-ectodomain          Cleaved, Asymmetric HIV-1 AE2 Env Ectodomain Trimer
         11      7txp  hiv-sosip-env-ectodomain3       Open, Symmetric D9/CD4-liganded HIV-1 SOSIP Env Ectodomain Trimer (ligands removed)
         12      5vn3  hiv-sosip-env-ectodomain4       Open, Symmetric 17b/CD4-liganded HIV-1 B41 SOSIP Env Ectodomain Trimer (ligands removed)
         13      2ins  insulin-hexamer                 DES-PHE B1 Bovine Insulin Hexamer
         14      4zxb  insulin-receptor-ectodomain     Human Insulin Receptor Ectodomain IRαβ
         15      7xix  sars-cov2-S-BA2                 Fully Glycosylated, Closed SARS-CoV-2 Omicron BA.2 Variant Spike
         16      6e8w  hiv-mpertm3-membrane1           HIV-1 Env MPER-TM Trimer in a DMPC/DHPC Symmetric Bilayer
         17      6e8w  hiv-mpertm3-membrane2           HIV-1 Env MPER-TM Trimer in an Asymmetric, Model Viral Bilayer
         18      5fkw  ecoli-polymerase                E. coli Replicative DNA Polymerase Complex Bound to DNA
         19      1mob  sperm-whale-myoglobin           Sperm whale myoglobin
         20    P22033  methylmalonyl-coa-mutase        Mitochondrial methylmalonyl-CoA mutase (Alphafold P22033)
         21      1aon  groel-groes-adp                 Asymmetric GroEL/GroES Chaperonin Complex
         22      2bgj  ferredoxin-fad                  Ferredoxin-NADP(H) Reductase from Rhodobacter capsulatus

.. _sub_fetch_example:

fetch-example
-------------

This subcommand is like ``run-example``, except it only copies the YAML input file needed to run the example to the CWD.  It can be edited, if desired, and then run using the ``run`` subcommand. 

