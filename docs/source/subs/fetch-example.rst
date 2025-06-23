fetch-example
-------------

Like ``run-example``, except it only copies the YAML input file needed to run the example to the CWD.  It can be edited, if desired, and then run using the ``run`` subcommand (see above).  If you request help on ``fetch-example``, you will see a list of examples:

.. code-block:: bash

   $ pestifer fetch-example --help
   usage: pestifer fetch-example [-h] [--config-updates CONFIG_UPDATES [CONFIG_UPDATES ...]] number

   Fetch YAML config for one of the examples:
   1: BPTI
   2: BPTI with phosphate ion excluded in salty solution
   3: BPTI, no phosphate, some random mutations plus deletion of one disulfide
   4: BPTI, no phosphate, introducing a disulfide via mutations
   5: HIV-1 protease dimer
   6: HIV-1 Env Trimer 4zmj
   7: HIV-1 Env Trimer 4tvp, Fabs and sulfate ions removed
   8: HIV-1 Env Trimer 8fad, drug molecule removed
   9: HIV-1 Env Trimer 8fae, drug molecule removed
   10: HIV-1 Env Trimer 7txd, substitutions, no ligands
   11: HIV-1 Env Trimer 5vn3, excluding sCD4 and Fab chains, with Gly3 stubs replacing missing v1/2
   12: hexameric insulin
   13: insulin receptor ectodomain, Fabs removed
   14: BA.2 SARS-CoV-2 Spike 7xix, fully glycosylated using grafts, and cleaved
   15: HIV-1 gp41 MPER-TM trimer 6e8w embedded in DMPC/DHPC bilayer
   16: HIV-1 gp41 MPER-TM trimer 6e8w embedded in model viral bilayer
   positional arguments:
   number                example number

   options:
   -h, --help            show this help message and exit
   --config-updates CONFIG_UPDATES [CONFIG_UPDATES ...]
                           yaml files to update example
