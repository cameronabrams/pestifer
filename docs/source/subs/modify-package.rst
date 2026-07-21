.. _subs modify-package:

modify-package
---------------

The subcommand is meant for development use only, and is only available when pestifer is installed as an **editable source package**. It allows modifications to your source package, including adding and deleting examples and updating atomselect macros.

This subcommand will only work on a full source repository, so if you want to use it, you will need to **fork the repository cameronabrams/pestifer on GitHub** and then clone your fork to your local machine.  If you have not already done so, you can do this with the following commands:

.. code-block:: bash

    git clone git@github.com:your-name/pestifer.git # assuming you have SSH access to your github account and you forked pestifer
    cd pestifer
    pip install -e . # so that you can tell `pestifer` to modify itself

.. _modify branching:

The Git workflow is automatic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because every modification is meant to become a pull request, ``modify-package`` folds the Git workflow into each command and **runs it by default**: it verifies your working tree is clean, makes the change on a **fresh branch** (auto-named ``modpkg/<category>-<verb>-<detail>``), and commits **exactly** the files it touched — then prints the ``git push`` / ``gh pr create`` steps.  You do not need to create a branch or commit by hand.

- ``--branch NAME`` — name the branch yourself instead of using the auto-generated name.
- ``--no-branch`` — skip branching and committing entirely; just apply the change to your working tree (useful when you are already on a feature branch, or want to stage the change yourself).

Both options are available on every ``example``, ``pdb-repo``, and ``charmmff`` verb.  The default (auto-branch) requires a clean working tree, so the branch contains only your contribution.

.. _modify add example:

Adding an Example
~~~~~~~~~~~~~~~~~

I developed most of the examples by iteration, so I would start with a simple YAML config file and then modify it to add new features or test new functionality.  The simplest way to start with a new example is to use ``pestifer new-system``.  For example, consider the :ref:`myoglobin example <example sperm-whale-myoglobin>`.  I started with:

.. code-block:: bash

    $ cd
    $ mkdir ~/1mob_working_directory # create a working directory
    $ cd ~/1mob_working_directory
    $ pestifer new-system 1mob # this is the PDB ID for sperm whale myoglobin
    $ pestifer build 1mob.yaml # this will run the psfgen task and generate a PSF file for the system

This run built successfully.  So I rebuilt a full template config:

.. code-block:: bash

    $ pestifer new-system 1mob --full

Then I edited ``1mob.yaml`` to add a salt concentration specification to the ``solvate`` task.  Then I ran it again:

.. code-block:: bash

    $ pestifer build 1mob.yaml

That also built successfully, so I added a new example to the pestifer package (from a clean working tree, so the auto-created branch contains only the new example):

.. code-block:: bash

    $ cd ~/1mob_working_directory  # where I have the 1mob.yaml file
    $ pestifer modify-package example add 1mob.yaml  # add the example to the package

This copies ``1mob.yaml`` to the appropriate location in the pestifer source package, updates the generated ``docs/source/examples/19/1mob.rst`` to link the new example, **and** — because the Git workflow is automatic (see :ref:`above <modify branching>`) — creates a branch ``modpkg/example-add-1mob-yaml`` and commits exactly those files.  Then I edited ``docs/source/examples/19/1mob.rst`` to add a description of the example, and folded that edit into the same commit:

.. code-block:: bash

    $ cd ~/Git/pestifer
    $ git add docs/source/examples/19/1mob.rst
    $ git commit --amend --no-edit

Then, I recompiled the documenation while my package was still on the ``modpkg/example-add-1mob-yaml`` branch:

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
    $ git merge modpkg/example-add-1mob-yaml
    $ git branch -d modpkg/example-add-1mob-yaml

If you want to add an example, you can do so in your own fork of the repository, and then submit a pull request to have it merged into the main repository.

.. _modify add residue:

Contributing a New Custom Residue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pestifer ships a small set of **built-in custom** residue definitions (for example, the CGenFF ligands ``83G`` and ``LF0``) that live in pestifer's shared ``charmmff/custom/`` directory (release-independent, a sibling of the per-release version directories) and supplement the standard CHARMM release.  You can contribute a new one from any file containing a CHARMM ``RESI`` block (a ``.str``, ``.rtf``, or ``.top`` file), and ``modify-package`` will:

1. validate the file and extract the ``RESI`` names it defines,
2. refuse names that already exist in the force field (unless you pass ``--force``),
3. copy the file into the shared ``charmmff/custom/`` directory,
4. register each ``RESI`` name under a segtype in :mod:`pestifer.core.labels` (default ``ligand``; use ``--segtype`` to choose another), and
5. clear the resource cache so the new residue is picked up on the next run.

The simplest invocation makes the change in your working tree:

.. code-block:: bash

    $ pestifer modify-package charmmff add-residue mylig.str --segtype ligand

Folding in the git workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Because a contribution is meant to become a pull request, the branch-and-commit step runs automatically (see :ref:`above <modify branching>`): from a clean working tree it creates a fresh branch off ``HEAD``, makes the change, and commits **exactly** the files it touched (the copied custom file and, if a segtype was registered, ``labels.py``) — never anything else in your tree.  Pass ``--branch NAME`` to name the branch, or ``--no-branch`` to apply the change without committing:

.. code-block:: bash

    $ pestifer modify-package charmmff add-residue mylig.str --branch add-mylig-residue

    Installed custom residue file: .../charmmff/custom/mylig.str
      RESI defined: MYLIG
      classified MYLIG as segtype 'ligand'
      resource cache cleared; it will rebuild on the next run.

    Committed to new branch 'add-mylig-residue':
        pestifer/resources/charmmff/custom/mylig.str
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

By default the collection is the residue's segtype (``lipid`` -> the ``lipid`` collection; ``ion``/``water`` -> ``solvent``); use ``--collection NAME`` to choose another, and ``--force`` to replace an entry already present for that resname.  As with the other verbs, the change is committed on a fresh branch by default (``--branch NAME`` to name it, ``--no-branch`` to skip) — from a clean working tree, containing **exactly** the changed collection tarball for you to push and open as a pull request.

Installing many entries at once
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``make-pdb-collection --streamID <stream>`` (or ``--substreamID``) writes a whole *directory* of entries -- one ``<RESI>/`` subdirectory per residue in the stream -- rather than a single entry.  ``pdb-repo add-entry`` accepts that container directory directly and installs every entry it finds in one command, instead of you having to invoke it once per residue:

.. code-block:: bash

    $ pestifer make-pdb-collection --streamID lipid --substreamID yeast --output-dir lipid-yeast
    # -> lipid-yeast/{DYPC/, DYPE/, PYPE/, YOPA/, ...}
    $ pestifer modify-package pdb-repo add-entry lipid-yeast --branch add-yeast-lipids

Each entry goes to its own default collection (its residue's segtype), so a mixed directory can populate several collections at once; ``--collection NAME`` forces them all into one.  Every affected tarball is validated up front and repacked once.  The argument is the same either way -- point ``add-entry`` at a single ``<RESI>/`` directory for one entry, or at a directory of them for a batch.

.. _modify regenerate segtypes:

Regenerating the derived segtype classification
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pestifer classifies most residues into a *segtype* (protein, lipid, glycan, nucleicacid, ligand, ion, water) by the CHARMM topology/stream file each residue is defined in, rather than from a hand-maintained list.  This classification is stored in the generated resource ``pestifer/resources/labels/derived_segtypes.json`` and merged into :attr:`Labels.segtype_of_resname <pestifer.core.labels.LabelMappers.segtype_of_resname>` at import.  After updating the bundled CHARMM force field (or adding a residue whose classification should follow from its defining file), regenerate it with

.. code-block:: bash

    $ pestifer modify-package charmmff regenerate-segtypes

This re-derives the classification from every installed CHARMM release, excluding the curated names in :mod:`pestifer.core.labels` (which always win), and rewrites the JSON.  Like the other developer actions, it can be combined with ``--branch`` to make the change on a fresh branch and commit it for review.


.. _modify ledger:

The modification ledger
~~~~~~~~~~~~~~~~~~~~~~~~~

Every mutating ``modify-package`` command records what it did in an append-only **ledger** at ``pestifer/resources/modifications.jsonl`` (one JSON object per line).  Each entry captures the category and verb, a one-line summary, the files touched, the committer, the branch, and a unique ``id``.  The ledger complements the resource *provenance* tags shown by :ref:`show-resources <subs_show_resources>` — provenance says a file is ``custom``; the ledger says *who added it, when, via which operation*.  Because the entry is committed alongside the change, it travels in the contribution's branch and PR.

List what has been recorded with ``ledger show``:

.. code-block:: bash

    $ pestifer modify-package ledger show            # all entries, oldest first
    $ pestifer modify-package ledger show --limit 10 --category charmmff

    #1  2026-07-08  charmmff/add-residue  add MYLIG to built-in custom [segtype: ligand]
    #2  2026-07-08  pdb-repo/add-entry    add MYLIG coordinates to the lipid collection

Reverting a modification
^^^^^^^^^^^^^^^^^^^^^^^^^

``ledger revert <id>`` reverses a recorded modification.  It finds the commit that made the change, ``git``-reverts it into the working tree, and then curates the ledger so the audit trail survives: the original entry is marked ``reverted`` and a new ``revert`` entry is appended.  Like the other verbs it runs on a fresh branch by default (``--branch NAME`` / ``--no-branch`` apply):

.. code-block:: bash

    $ pestifer modify-package ledger revert 1

    Reverted modification #1 (add MYLIG to built-in custom [segtype: ligand]); recorded as #3.

    Committed the revert to new branch 'modpkg/revert-1':
    ...

Revert uses ``git`` under the hood, so it works uniformly for any verb and surfaces a conflict if a later change touched the same files.  Only modifications that were **committed** (the default flow) can be reverted this way; a change made with ``--no-branch`` is left to you to undo by hand.
