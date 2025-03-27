.. _subs_run:

run
---

Builds a system according to the input configuration file.

.. code-block:: console

   $ pestifer run <config.yaml>

Here ``config.yaml`` is the name of the configuration file that describes the build.  Minimally, a pestifer config file for a run must have a ``tasks`` directive that specifies the *ordered list* of tasks the run should perform.  Typically, each tasks inherits a PSF/PDB/COOR/XSC NAMD file set for a system from a previous task.  

Detailed explanation of some *selected* common tasks you can use is below.

.. toctree:: 
   :maxdepth: 1

   runtasks/psfgen
   runtasks/md
   runtasks/ligate
   runtasks/solvate
   runtasks/bilayer
   runtasks/terminate
   runtasks/cleave

Please consult the Config Reference pages :ref:`config_ref tasks` for a full list of available pestifer run tasks.

We provide several :ref:`examples` that show a variety of task lists.
