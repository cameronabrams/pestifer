.. _subs_build:

build
-----

The main purpose of Pestifer is to prepare simulation-ready systems.  The ``build`` subcommand is used to run a system preparation from a configuration file.  The configuration file describes the tasks that should be performed to prepare the system, and the order in which they should be performed.

.. note::

   ``run`` is accepted as a synonym for ``build`` for backwards compatibility.

.. code-block:: bash

   $ pestifer build config.yaml

Here ``config.yaml`` is the name of the configuration file that describes the system preparation.  Minimally, a pestifer config file for a build must have a ``tasks`` directive that specifies the *ordered list* of tasks the build should perform.  

The first task must be one of a :ref:`fetch <subs_buildtasks_fetch>` or a :ref:`continuation <subs_buildtasks_continuation>` or a :ref:`merge <subs_buildtasks_merge>`.  ``fetch`` is used when you begin a run from an externally provided structure file (a PDB or mmCIF/PDBx-format file), while ``continuation`` is used when you begin a run from an already prepared system state, represented by PDB, PSF, and XSC files. ``merge`` is used when you want to combine two or more system states into a single state, for example to combine a protein structure with a ligand structure.  The final task must be a :ref:`terminate <subs_buildtasks_terminate>`, which is responsible for writing out the final system state and any associated files, and optionally packaging the system for production elsewhere.  In between, you can have any number of tasks that perform various operations on the system state, such as building the PSF, running minimization or MD simulations, or performing various manipulations of the system (e.g. ligation, cleavage, solvation, etc.).

Generally, each ``task`` is responsible for one or more changes of the system **state**, and Pestifer ensures that every task knows the current state of the system when it begins.  Here is an example of a task list that begins from the X-ray crystal structure of the bovine pancreatic trypsin inhibitor structure (PDB ID 6pti), builds the PSF and PDB, then subjects that state
to minimization and a short NVT MD simulation, and finally packaging for production elsewhere.

.. code-block:: yaml

   tasks:
     - fetch:
        - sourceID: 6pti
     - psfgen:
     - md:
         ensemble: minimize
     - md:
         ensemble: NVT
     - terminate:
         basename: my_system
         package:
           basename: prod_system
           md:
             ensemble: NPT

``pestifer build`` begins by parsing the configuration file, and then it sets up a ``Controller`` object that manages the tasks.  Before running anything, the controller statically **validates the task list**: each task declares the pipeline "currencies" it requires and provides (a fetched coordinate source, the built ``state`` fileset, the in-memory molecule object, molecular-dynamics output), and the controller checks the hand-offs so malformed pipelines fail immediately with an explanation rather than partway through a long build.  It flags, for example, a transform task with no system to act on, a ``psfgen`` that builds from a fetched source with no preceding ``fetch``, a task that needs the in-memory molecule object (``cleave``/``ligate``/``pdb2pqr``) after a ``continuation`` or ``merge`` (which provide files, but not that object), an ``mdplot`` with no preceding ``md``, or a task placed after the terminal ``terminate``.  The controller then executes each task in the order specified in the configuration file; if a task fails, the controller stops the build and reports the error.

Because the contract is pipeline-aware, ``psfgen`` **infers its mode from what precedes it**: after a ``fetch`` it builds a new topology from the fetched coordinates, while after a ``continuation`` or ``merge`` that established a state with no fetched source, it *edits the incoming topology* instead of rebuilding it.  See :ref:`psfgen <subs_buildtasks_psfgen>` for what that mode can and cannot do.

.. mermaid:: 
   :caption: Execution flow for the pestifer controller

   graph TD;
      A[Config] --> C{Next task?};
      subgraph group Controller [Controller];
         C --> Yes --> D[Run task];
         D -->|Success| C;
         D -->|Failure| F[Report error];
      end
      C --> No --> E[End];
      F --> E;

Checking a config before building
=================================

Everything the controller validates before it runs anything -- the config parses, every key is
one the schema defines, ``vmd``/``namd3``/``charmrun``/``catdcd`` resolve on your PATH and are
executable, the task hand-offs are sound -- can be checked without starting a build:

.. code-block:: console

  $ pestifer build config.yaml --check

This reports the resolved toolchain, the CHARMM force-field release, the NAMD mode and
processor count, and the task plan the controller would execute, then exits.  It runs in about
a second, fetches nothing, and writes nothing:

.. code-block:: text

   pestifer 3.16.2 build --check: 12er.yaml

     title       New template pestifer config for id 12er (PDB)
     charmmff    feb26
     namd        cpu mode, 24 PEs
     charmrun    /usr/local/bin/charmrun
     namd3       /usr/local/bin/namd3
     vmd         /usr/local/bin/vmd
     catdcd      /usr/local/bin/catdcd

     task plan (8 tasks):
        0  fetch
        1  psfgen
        2  md
        3  solvate
        4  md
        5  md
        6  density_equilibrate
        7  terminate

   check passed; this config would build

The exit status is 0 when the build would start and 1 when it would not, so ``--check`` can gate
a batch script.  A mistyped directive is reported with the alternatives the schema does accept:

.. code-block:: text

   ERROR: Attribute 'mutationz' invalid; expecting one of ['mutations', 'ssbonds',
   'ssbondsdelete', 'links', 'deletions', 'substitutions', ...] under 'mods'.

   check FAILED; the build would not start

Two conditions are reported as **warnings** rather than failures: a config value that names an
input file not present in the working directory, and a working directory that already contains
other files (a build writes a great many, and expects a directory of its own).  Neither stops
the build.

Adding ``--json`` writes the same report to standard output as a JSON object -- the resolved
executables, the task list, and the ``warnings`` and ``errors`` arrays -- while the ordinary log
output stays on standard error.  This is the form to use when a script or an agent is driving
pestifer; see :ref:`agent_driven_builds`.

Resuming an interrupted build
=============================

Long builds -- a membrane system with hours of relaxation MD, say -- can be interrupted by a wall-clock limit, a node failure, or a Ctrl-C.  As it runs, pestifer writes a **run manifest** (``.pestifer-manifest.json``) into the run directory recording, for each task, the spec it was given and the state it produced.  You can resume from it:

.. code-block:: console

  $ pestifer build config.yaml --restart

``--restart`` reads the manifest and resumes from the last **cleanly-completed** task, discarding any partial files the interrupted task left behind.  Tasks are matched by a hash of their specs taken *before* execution, so if you edited the config between runs, the first task whose spec changed becomes the resume point and everything after it re-runs.

.. code-block:: console

  $ pestifer build config.yaml --from make_membrane_system   # resume from a task you name
  $ pestifer build config.yaml --fresh                       # ignore any manifest; build from scratch

``--from`` takes a task index or task name and overrides the auto-detected resume point (it implies ``--restart``).  ``--fresh`` does the opposite, ignoring an existing manifest entirely.

.. note::

   A resume point must be a clean hand-off of the pipeline ``state``.  Pestifer refuses to resume into a task whose predecessor did not produce a state fileset, rather than silently starting from an incomplete system.

Detailed explanation of some *selected* common tasks you can use is below.

.. toctree::
   :maxdepth: 1

   buildtasks/fetch
   buildtasks/continuation
   buildtasks/merge
   buildtasks/psfgen
   buildtasks/validate
   buildtasks/pdb2pqr
   md <buildtasks/mdtask>
   buildtasks/density_equilibrate
   buildtasks/membrane_equilibrate
   buildtasks/manipulate
   buildtasks/ligate
   buildtasks/cleave
   buildtasks/solvate
   buildtasks/make_membrane_system
   buildtasks/ring_check
   buildtasks/mdplot
   buildtasks/terminate

Please consult the Configuration Reference pages for :ref:`tasks <config_ref tasks>`, for a full list of available pestifer build tasks.

We provide several :ref:`examples` that show a variety of task lists.
