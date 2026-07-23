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

``pestifer build`` begins by parsing the configuration file, and then it sets up a ``Controller`` object that manages the tasks.  Before running anything, the controller statically **validates the task list**: each task declares the pipeline "currencies" it requires and provides (a fetched coordinate source, the built ``state`` fileset, the in-memory molecule object, molecular-dynamics output), and the controller checks the hand-offs so malformed pipelines fail immediately with an explanation rather than partway through a long build.  It flags, for example, a transform task with no system to act on, a ``continuation`` followed by ``psfgen`` (which would rebuild from scratch and discard the continued system), an ``mdplot`` with no preceding ``md``, or a task placed after the terminal ``terminate``.  The controller then executes each task in the order specified in the configuration file; if a task fails, the controller stops the build and reports the error.

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
   buildtasks/manipulate
   buildtasks/ligate
   buildtasks/cleave
   buildtasks/solvate
   buildtasks/make_membrane_system
   buildtasks/ring_check
   buildtasks/terminate

Please consult the Configuration Reference pages for :ref:`tasks <config_ref tasks>`, for a full list of available pestifer build tasks.

We provide several :ref:`examples` that show a variety of task lists.
