Usage
=====

Installation of the ``pestifer`` package gives access to the ``pestifer`` command.  The general syntax for invoking ``pestifer`` is

.. code-block:: console

   $ pestifer <subcommand> <options>

Help with any subcommand can be obtained via

.. code-block:: console

   $ pestifer <subcommand> --help

General options for all subcommands:

.. code-block:: console

  --no-banner          turn off the banner
  --loglevel LOGLEVEL  Log level for messages written to diagnostic log (debug|info)
  --diag DIAG          diagnostic log file

There are several subcommands:

``run``: the main subcommand that uses a user's input config file to build a system.

.. code-block:: console

   $ pestifer run <config.yaml>

Here ``config.yaml`` is the name of the configuration file that describes the build.  The
best way to learn about the configuration file is to run the examples provided, and
then you can try to make your own.

``run-example``: there are 21 example systems; to run number four, for example:

.. code-block:: console
   
   $ pestifer run-example 4

(Best to do that in a clean directory.)

``config-help``: Interactive help in constructing a config file.

.. code-block:: console

   $ pestifer config-help
   Help on user-provided configuration file format
       Help available for charmmff, psfgen, namd2, title, paths, tasks
   $ pestifer config-help tasks
   Help on user-provided configuration file format
   tasks:
      Specifies the tasks to be performed; each is a dictionary with a
         heading which is a reserved task name
      type: list
      Help available for psfgen, ligate, solvate, relax, terminate
   $ pestifer config-help tasks psfgen
   Help on user-provided configuration file format
   tasks->
   psfgen:
      Parameters controlling initial psfgen run
      type: dict
      Help available for source, mods, minimize, cleanup

