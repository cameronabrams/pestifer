.. _subs_run:

run
---

The main purpose of Pestifer is to prepare simulation-ready systems.  The ``run`` subcommand is used to run a system preparation from a configuration file.  The configuration file describes the tasks that should be performed to prepare the system, and the order in which they should be performed.

.. code-block:: bash

   $ pestifer run config.yaml

Here ``config.yaml`` is the name of the configuration file that describes the system preparation.  Minimally, a pestifer config file for a run must have a ``tasks`` directive that specifies the *ordered list* of tasks the run should perform.  

The first task must either be a :ref:`fetch <subs_runtasks_fetch>` or a :ref:`continuation <subs_runtasks_continuation>`.  ``fetch`` is used when you begin a run from an externally provided structure file (a PDB or mmCIF/PDBx-format file), while ``continuation`` is used when you begin a run from an already prepared system state, represented by PDB, PSF, and XSC files.

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

``pestifer run`` begins by parsing the configuration file, and then it sets up a ``Controller`` object that manages the tasks.  The controller executes each task in the order specified in the configuration file, and it handles any dependencies between tasks.  If a task fails, the controller will stop the run and report the error.

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

   runtasks/fetch
   runtasks/continuation
   runtasks/psfgen
   runtasks/validate
   runtasks/pdb2pqr
   md <runtasks/mdtask>
   runtasks/ligate
   runtasks/cleave
   runtasks/solvate
   runtasks/make_membrane_system
   runtasks/terminate

Please consult the Configuration Reference pages for :ref:`tasks <config_ref tasks>`, for a full list of available pestifer run tasks.

We provide several :ref:`examples` that show a variety of task lists.
