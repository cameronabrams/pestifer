.. _subs modify-package:

modify-package
---------------

The subcommand is meant for development use only. It allows modifications to the **source package**, including adding and deleting examples and updating atomselect macros.

This subcommand will only work on a full source repository, so if you want to use it, you will need to fork the repository on GitHub and then clone your fork to your local machine.  If you have not already done so, you can do this with the following commands:

.. code-block:: bash

    git clone git@github.com:your-name/pestifer.git # assuming you have SSH access to your github account and you forked pestifer
    cd pestifer
    pip install -e . # so that you can tell `pestifer` to modify itself

I also recommend creating a Git branch for your modifications, so you can easily revert them if needed:

.. code-block:: bash

    git checkout -b my-modifications

Adding an Example
~~~~~~~~~~~~~~~~~

I developed most of the examples by iteration, so I would start with a simple YAML config file and then modify it to add new features or test new functionality.  The simplest way to start with a new example is to use ``pestifer new-system``.  For example, consider the :ref:`myoglobin example <example 1mob>`.  I started with:

.. code-block:: bash

    $ cd
    $ mkdir ~/1mob_working_directory # create a working directory
    $ cd ~/1mob_working_directory
    $ pestifer new-system --id 1mob # this is the PDB ID for sperm whale myoglobin
    $ pestifer run 1mob.yaml # this will run the psfgen task and generate a PSF file for the system

This run built successfully.  So I rebuilt a full template config:

.. code-block:: bash

    $ pestifer run 1mob.yaml --full

Then I edited ``1mob.yaml`` to add a salt concentration specification to the ``solvate`` task.  Then I ran it again:

.. code-block:: bash

    $ pestifer run 1mob.yaml

That also built successfully, so I added a new example to the ``examples`` directory:

.. code-block:: bash

    $ cd ~/Git/pestifer        # my local source repository
    $ git branch -b new-1mob   # create a new branch for the example
    $ cd ~/1mob_working_directory  # where I have the 1mob.yaml file
    $ pestifer modify-package --add-example 1mob.yaml  # add the example to the package

This will copy the file ``1mob.yaml`` to the ``examples`` directory in the pestifer source package, and it will also update the ``docs/source/examples/1mob.rst`` file to include a link to the new example.  Then I edited the ``docs/source/examples/1mob.rst`` file to add a description of the example and how it works.  Finally, I committed the changes:

.. code-block:: bash

    $ git commit -m "Add 1mob example"

Then, I recompiled the documenation:

.. code-block:: bash

    $ cd ~/Git/pestifer/docs
    $ make html

And I made sure the new example was there.  It also appeared correctly as number 19 when displaying the examples:

.. code-block:: bash

    $ pestifer show-resources --examples

        Examples:

        Index        ID  Name                            Description
            1      6pti  bpti1                           Bovine Pancreatic Trypsin Inhibitor (BPTI)
            2      6pti  bpti2                           BPTI Excluding the Phosphate Ion
            3      6pti  bpti3                           BPTI with a Mutated-out Disulfide Bond
            4      6pti  bpti4                           BPTI with a Mutated-in Disulfide Bond
            5      1f7a  hiv-protease                    HIV Protease with Patches to Protonate Aspartates
            6      1fas  toxin                           Green Mamba Toxin at pH 7.0
            7      4zmj  hiv-sosip-env-ectodomain1       Closed, Unliganded HIV-1 BG505 Env SOSIP.664 Trimer
            8      4tvp  hiv-sosip-env-ectodomain2       Closed, PGT122/35O22-Liganded HIV-1 BG505 Env SOSIP.664 Trimer (ligands removed)
            9      8fad  hiv-ad8-env-ectodomain          Cleaved, Asymmetric HIV-1 AD8 Env Ectodomain Trimer
           10      8fae  hiv-ae2-env-ectodomain          Cleaved, Asymmetric HIV-1 AE2 Env Ectodomain Trimer
           11      7txd  hiv-sosip-env-ectodomain3       Open, Symmetric D9/CD4-liganded HIV-1 SOSIP Env Ectodomain Trimer (ligands removed)
           12      5vn3  hiv-sosip-env-ectodomain4       Open, Symmetric 17b/CD4-liganded HIV-1 B41 SOSIP Env Ectodomain Trimer (ligands removed)
           13      2ins  insulin-hexamer                 DES-PHE B1 Bovine Insulin Hexamer
           14      4zxb  insulin-receptor-ectodomain     Human Insulin Receptor Ectodomain IRαβ
           15      7xix  sars-cov2-S-BA2                 Fully Glycosylated, Closed SARS-CoV-2 Omicron BA.2 Variant Spike
           16      6e8w  hiv-mpertm3-membrane1           HIV-1 Env MPER-TM Trimer in a DMPC/DHPC Symmetric Bilayer
           17      6e8w  hiv-mpertm3-membrane2           HIV-1 Env MPER-TM Trimer in an Asymmetric, Model Viral Bilayer
           18      5fkw  ecoli-polymerase                E. coli Replicative DNA Polymerase Complex Bound to DNA
           19      1mob  1mob                            Sperm whale myoglobin

Satisfied with the example, I merged the branch back into ``main``:

.. code-block:: bash

    $ git checkout main
    $ git merge new-1mob
    $ git branch -d new-1mob

If you want to add an example, you can do so in your own fork of the repository, and then submit a pull request to have it merged into the main repository.

