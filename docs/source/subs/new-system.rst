.. _subs new-system:

new-system
-----------------

Pestifer provides a subcommand ``new-system`` as a convenient way to generate a bare-bones system configuration file for a new simulation.  This is useful when you want to start a new simulation from scratch.

``pestifer new-system`` accepts a single argument, which is intepreted as either the PDB ID or AlphaFold ID of the new system.

For example, to create a new system configuration file for the PDB ID 1abc, you would run:

.. code-block:: bash

   $ pestifer new-system 1abc

This will create a new file named ``1abc.yaml`` in the current directory, containing a basic configuration for the system that, by default, has only a ``fetch`` task and a single ``psfgen`` task:

.. code-block:: yaml

    title: New template pestifer input configuration for PDB ID 1abc
    tasks:
      - fetch:
          sourceID: 1abc
      - psfgen:
           
Using this config is helpful in getting over the first hurdle of a new system, which is making sure pestifer can write a good, albeit minimal, ``psfgen`` script.

The ``new-system`` subcommand bases its output on the example file for example 1, the simple :ref:`BPTI system <example bpti1>` .

There are several options to ``new-system``:

``--full``
+++++++++++

Including this flag option will add solvation, minimization, and equilibration tasks to the configuration file.  For example:

.. code-block:: bash

   $ pestifer new-system 1abc --full

This will create a new file named ``1abc.yaml`` in the current directory, containing a more complete configuration for the system, including tasks for solvation, minimization, and equilibration.

``--output <filename>``
++++++++++++++++++++++++

This option allows you to specify the name of the output file for the new system configuration.  For example:

.. code-block:: bash

   $ pestifer new-system 1abc --output my_custom_config.yaml

This will create a new file named ``my_custom_config.yaml`` in the current directory, containing the same configuration as ``1abc.yaml``, but with the specified filename.

``--title <title>``
+++++++++++++++++++++++

This option allows you to specify a custom title for the new system configuration.  For example:

.. code-block:: bash

   $ pestifer new-system 1abc --title "My Custom Title"

This will create a new file named ``1abc.yaml`` in the current directory with the specified title.
