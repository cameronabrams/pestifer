.. _usage:

Usage
=====

Pestifer's main functionality is available on the command line.  You can also access Pestifer's library of useful TcL/VMD procs from your own VMD scripts.

Command-line Usage
------------------

Installation of Pestifer gives access to the ``pestifer`` command.  The general syntax for invoking ``pestifer`` is

.. code-block:: bash

   $ pestifer <subcommand> <options>

Help with any subcommand can be obtained via

.. code-block:: bash

   $ pestifer <subcommand> --help

General options for all subcommands:

.. code-block:: bash

  --no-banner          turn off the banner
  --loglevel LOGLEVEL  Log level for messages written to diagnostic log (debug|info)
  --diag DIAG          diagnostic log file

Subcommands
-----------

.. toctree::
   :maxdepth: 1

   subs/config-help
   subs/build
   subs/build-example
   subs/new-system
   subs/wheretcl
   subs/make-pdb-collection
   subs/desolvate
   subs/make-namd-restart
   subs/show-resources
   subs/follow-namd-log
   subs/mdplot
   subs/modify-package

.. _use in vmd scripts:

Use in VMD scripts
------------------

Pestifer provides a library of useful Tcl/VMD procs and VMD macro definitions.
To make them available in every VMD session, run once (from the conda environment
in which pestifer is installed):

.. code-block:: bash

   $ pestifer setup-vmd

This writes ``~/.pestifer/vmd_init.tcl`` with the absolute path to pestifer's Tcl
library hardcoded, and appends a one-line source hook to ``~/.vmdrc``:

.. code-block:: tcl

   # pestifer setup-vmd
   if {[file exists ~/.pestifer/vmd_init.tcl]} { source ~/.pestifer/vmd_init.tcl }

That single line is the only thing that ever needs to be in ``~/.vmdrc`` — it does
not change between pestifer versions.  After upgrading pestifer or switching conda
environments, simply re-run ``pestifer setup-vmd`` to regenerate
``~/.pestifer/vmd_init.tcl`` with the updated path.

.. important::

   Run ``pestifer setup-vmd`` from the shell environment (or conda environment)
   in which pestifer is installed, so that the correct installation path is
   recorded.  VMD itself can be launched from any environment — no dynamic
   ``exec pestifer`` call is made at VMD startup.

Once set up, VMD automatically loads pestifer's Tcl library on startup, including
extended definitions of the ``glycan`` and ``lipid`` atomselect macros and all
PestiferUtil procedures.