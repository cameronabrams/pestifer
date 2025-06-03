make-namd-restart
-----------------

Pestifer provides a subcommand ``make-namd-restart`` as a convenient way to generate inputs for sequential NAMD runs, or to complete NAMD runs that terminated early.

For example, say you have a run whose config file was ``my_run_01.namd`` and which generated successful output files ``my_run_01.coor``, ``my_run_01.vel``, ``my_run_01.xsc``, etc.  Suppose the log file is ``my_run_01.log``.  In such a case, we refer to the string ``my_run_01`` as the "basename" of this run; it is also the argument of ``outputname`` in a NAMD config.  ``pestifer make-namd-restart`` can be used to generate the next config file like this:

.. code-block:: console

    $ pestifer make-namd-restart --log my_run_01.log --config my_run_01.namd --new-base my_run_02 --run 1000000

This will create the new config file ``my_run_02.namd`` and it will specify a run of 1,000,000 timesteps, with inputs comprising the successful outputs defined in ``my_run_01.namd``.

If you are running in a SLURM batch environment, and your SLURM ``bash`` script contains an assignment to the variable ``BASENAME``, you can include that as an option to ``pestifer make-namd-restart`` and it will update the SLURM script for you; e.g.: 

.. code-block:: console

    $ pestifer make-namd-restart --log my_run_01.log --config my_run_01.namd --new-base my_run_02 --run 1000000 --slurm my_slurm.sh

This will replace an assignment statement to the variable ``BASENAME`` with

.. code-block:: bash

    BASENAME=my_run_02
