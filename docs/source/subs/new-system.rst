.. _subs new-system:

new-system
-----------------

Pestifer provides a subcommand ``new-system`` as a convenient way to generate a bare-bones system configuration file for a new simulation.  This is useful when you want to start a new simulation from scratch, or when you want to create a new system based on an existing one.

``pestifer new-system`` accepts a single argument, which is intepreted as either the PDB ID or AlphaFold ID of the new system.

For example, to create a new system configuration file for the PDB ID 1abc, you would run:

.. code-block:: bash

   $ pestifer new-system --id 1abc

This will create a new file named `1abc.yaml` in the current directory, containing a basic configuration for the system that, by default,
has only a single ``psfgen`` task:

.. code-block:: yaml

    title: New template pestifer input configuration for PDB ID 1abc
    tasks:
      - psfgen:
          id: 1abc
           
Using this config is helpful in getting over the first hurdle of a new system, which is making sure pestifer can write a good psfgen script.

The ``--full`` option will add solvation, minimization, and equilibration tasks to the configuration file.  For example:

.. code-block:: bash

   $ pestifer new-system --id 1abc --full

This will create a new file named `1abc.yaml` in the current directory, containing a more complete configuration for the system, including tasks for solvation, minimization, and equilibration.