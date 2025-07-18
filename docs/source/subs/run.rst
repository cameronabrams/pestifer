.. _subs_run:

run
---

The main purpose of Pestifer is to prepare simulation-ready systems.  The ``run`` subcommand is used to run a system preparation from a configuration file.  The configuration file describes the tasks that should be performed to prepare the system, and the order in which they should be performed.

.. code-block:: bash

   $ pestifer run config.yaml

Here ``config.yaml`` is the name of the configuration file that describes the system preparation.  Minimally, a pestifer config file for a run must have a ``tasks`` directive that specifies the *ordered list* of tasks the run should perform.  Typically, each tasks inherits a PSF/PDB/COOR/XSC NAMD file set for a system from a previous task.  For example:

.. code-block:: yaml

   tasks:
     - psfgen:
         source:
           id: 6pti
     - md:
         ensemble: minimize
     - md:
         ensemble: NVT
     - terminate:
         basename: my_system
         package:
           ensemble: NPT
           basename: prod_system

``pestifer run`` begins by parsing the configuration file, and then it sets up a ``Controller`` object that manages the tasks.  The controller executes each task in the order specified in the configuration file, and it handles any dependencies between tasks.  If a task fails, the controller will stop the run and report the error.

.. mermaid:: 
   :caption: Pestifer Run Task Flowchart

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

   runtasks/psfgen
   runtasks/pdb2pqr
   runtasks/md
   runtasks/cleave
   runtasks/ligate
   runtasks/solvate
   runtasks/make_membrane_system
   runtasks/terminate

Please consult the Configuration Reference pages for the ``tasks`` directive, :ref:`config_ref tasks`, for a full list of available pestifer run tasks.

We provide several :ref:`examples` that show a variety of task lists.
