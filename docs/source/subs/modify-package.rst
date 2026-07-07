.. _subs modify-package:

modify-package
---------------

The subcommand is meant for development use only, and is only available when pestifer is installed as an **editable source package**. It allows modifications to your source package, including adding and deleting examples and updating atomselect macros.

This subcommand will only work on a full source repository, so if you want to use it, you will need to **fork the repository cameronabrams/pestifer on GitHub** and then clone your fork to your local machine.  If you have not already done so, you can do this with the following commands:

.. code-block:: bash

    git clone git@github.com:your-name/pestifer.git # assuming you have SSH access to your github account and you forked pestifer
    cd pestifer
    pip install -e . # so that you can tell `pestifer` to modify itself

I also recommend creating a Git branch for your modifications, so you can easily revert them if needed:

.. code-block:: bash

    git checkout -b my-modifications

.. _modify add example:

Adding an Example
~~~~~~~~~~~~~~~~~

I developed most of the examples by iteration, so I would start with a simple YAML config file and then modify it to add new features or test new functionality.  The simplest way to start with a new example is to use ``pestifer new-system``.  For example, consider the :ref:`myoglobin example <example sperm-whale-myoglobin>`.  I started with:

.. code-block:: bash

    $ cd
    $ mkdir ~/1mob_working_directory # create a working directory
    $ cd ~/1mob_working_directory
    $ pestifer new-system --id 1mob # this is the PDB ID for sperm whale myoglobin
    $ pestifer build 1mob.yaml # this will run the psfgen task and generate a PSF file for the system

This run built successfully.  So I rebuilt a full template config:

.. code-block:: bash

    $ pestifer build 1mob.yaml --full

Then I edited ``1mob.yaml`` to add a salt concentration specification to the ``solvate`` task.  Then I ran it again:

.. code-block:: bash

    $ pestifer build 1mob.yaml

That also built successfully, so I added a new example to the pestifer package **after** switching to a new branch:

.. code-block:: bash

    $ cd ~/Git/pestifer        # my local source repository
    $ git branch -b new-1mob   # create a new branch for the example
    $ cd ~/1mob_working_directory  # where I have the 1mob.yaml file
    $ pestifer modify-package example add 1mob.yaml  # add the example to the package

This will copy the file ``1mob.yaml`` to the appropriate location in pestifer source package, and it will also update the ``docs/source/examples/19/1mob.rst`` file to include a link to the new example.  Then I edited the ``docs/source/examples/19/1mob.rst`` file to add a description of the example and how it works.  Finally, I committed the changes:

.. code-block:: bash

    $ git add .
    $ git commit -m "Add 1mob example"

Then, I recompiled the documenation while my package was still in the ``new-1mob`` branch:

.. code-block:: bash

    $ cd ~/Git/pestifer/docs
    $ make html

And I made sure the new example was there.  It also appeared correctly as number 19 when displaying the examples:

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
           16      6e8w  hiv-mpertm3-membrane1           HIV-1 gp41 MPER-TM trimer 6e8w embedded in DMPC/DHPC bilayer
           17      6e8w  hiv-mpertm3-membrane2           HIV-1 gp41 MPER-TM trimer 6e8w embedded in model viral bilayer
           18      5fkw  ecoli-polymerase                E. Coli replicative DNA polymerase complex bound to a primer-template DNA
           19      1mob  1mob                            Sperm whale myoglobin

I then renamed this example using the ``modify-package`` subcommand:

.. code-block:: bash

    $ pestifer modify-package example rename 19 sperm-whale-myoglobin

Satisfied with the example, I merged the branch back into ``main``:

.. code-block:: bash

    $ git checkout main
    $ git merge new-1mob
    $ git branch -d new-1mob

If you want to add an example, you can do so in your own fork of the repository, and then submit a pull request to have it merged into the main repository.

.. _modify add residue:

Contributing a New Custom Residue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pestifer ships a small set of **built-in custom** residue definitions (for example, the CGenFF ligands ``83G`` and ``LF0``) that live in the active force field's ``custom/`` directory and supplement the standard CHARMM release.  You can contribute a new one from any file containing a CHARMM ``RESI`` block (a ``.str``, ``.rtf``, or ``.top`` file), and ``modify-package`` will:

1. validate the file and extract the ``RESI`` names it defines,
2. refuse names that already exist in the force field (unless you pass ``--force``),
3. copy the file into the force field's ``custom/`` directory,
4. register each ``RESI`` name under a segtype in :mod:`pestifer.core.labels` (default ``ligand``; use ``--segtype`` to choose another), and
5. clear the resource cache so the new residue is picked up on the next run.

The simplest invocation makes the change in your working tree:

.. code-block:: bash

    $ pestifer modify-package charmmff add-residue mylig.str --segtype ligand

Folding in the git workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because a contribution is meant to become a pull request, the ``--branch`` option folds the branch-and-commit step into the command.  It **requires a clean working tree** (so the resulting branch contains only your contribution), creates a new branch off your current ``HEAD``, makes the change, and commits **exactly** the files it touched (the copied custom file and, if a segtype was registered, ``labels.py``) — never anything else in your tree:

.. code-block:: bash

    $ pestifer modify-package charmmff add-residue mylig.str --branch add-mylig-residue

    Installed custom residue file: .../charmmff/feb26/custom/mylig.str
      RESI defined: MYLIG
      classified MYLIG as segtype 'ligand'
      resource cache cleared; it will rebuild on the next run.

    Committed to new branch 'add-mylig-residue':
        pestifer/resources/charmmff/feb26/custom/mylig.str
        pestifer/core/labels.py

    Review it with `git show`, then push and open a pull request:
        git push -u origin add-mylig-residue
        gh pr create      # or open the PR on GitHub

Pushing the branch and opening the pull request are left to you — inspect the commit with ``git show`` first, then push to your fork and open a PR against ``cameronabrams/pestifer`` for review.

.. _modify add pdb entry:

Contributing PDB-repository coordinates for a residue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Adding a residue's *topology* (:ref:`above <modify add residue>`) lets pestifer build the residue, but ``make_membrane_system`` also needs sampled *coordinates* to place a lipid on the grid.  Those live in pestifer's built-in **PDB repository**: a set of *collections*, each a ``<stream>.tgz`` tarball whose single top-level directory is the stream name and which holds one ``<RESI>/`` subdirectory per residue (its ``info.yaml`` plus a set of conformer PDBs).

First generate the coordinates with :ref:`make-pdb-collection <sub_make_pdb_collection>`, which samples conformers of the residue and writes the entry directory:

.. code-block:: bash

    $ pestifer make-pdb-collection --resname MYLIP --output-dir mycoords
    # -> mycoords/MYLIP/{info.yaml, MYLIP-00.pdb, MYLIP-01.pdb, ...}

Then install that entry into the repository with ``pdb-repo add-entry``, which validates it (``info.yaml`` parses and every conformer PDB it names is present) and adds it under ``<collection>/MYLIP/`` in the target collection tarball, creating the tarball if the collection does not exist yet:

.. code-block:: bash

    $ pestifer modify-package pdb-repo add-entry mycoords/MYLIP --branch add-mylip-coords

By default the collection is the residue's segtype (``lipid`` -> the ``lipid`` collection; ``ion``/``water`` -> ``water_ions``); use ``--collection NAME`` to choose another, and ``--force`` to replace an entry already present for that resname.  As with the other contribution flows, ``--branch`` requires a clean working tree and commits **exactly** the changed collection tarball on a new branch for you to push and open as a pull request.

.. _modify regenerate segtypes:

Regenerating the derived segtype classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pestifer classifies most residues into a *segtype* (protein, lipid, glycan, nucleicacid, ligand, ion, water) by the CHARMM topology/stream file each residue is defined in, rather than from a hand-maintained list.  This classification is stored in the generated resource ``pestifer/resources/labels/derived_segtypes.json`` and merged into :attr:`Labels.segtype_of_resname <pestifer.core.labels.LabelMappers.segtype_of_resname>` at import.  After updating the bundled CHARMM force field (or adding a residue whose classification should follow from its defining file), regenerate it with

.. code-block:: bash

    $ pestifer modify-package charmmff regenerate-segtypes

This re-derives the classification from every installed CHARMM release, excluding the curated names in :mod:`pestifer.core.labels` (which always win), and rewrites the JSON.  Like the other developer actions, it can be combined with ``--branch`` to make the change on a fresh branch and commit it for review.

