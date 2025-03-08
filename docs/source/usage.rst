.. _usage:

Usage
=====

Pestifer's main functionality is available on the command line.  You can also access pestifer's library of useful TcL/VMD procs from your own VMD scripts.

Command-line Usage
------------------

Installation of pestifer gives access to the ``pestifer`` command.  The general syntax for invoking ``pestifer`` is

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

Subcommands
-----------

.. toctree::
   :maxdepth: 1

   subs/config-help
   subs/run
   subs/run-example
   subs/fetch-example
   subs/wheretcl
   subs/inittcl
   subs/make-resi-database
   subs/desolvate
   subs/make-namd-restart
   subs/show-resources

