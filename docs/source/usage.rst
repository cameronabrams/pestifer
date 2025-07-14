.. _usage:

Usage
=====

Pestifer's main functionality is available on the command line.  You can also access pestifer's library of useful TcL/VMD procs from your own VMD scripts.

Command-line Usage
------------------

Installation of pestifer gives access to the ``pestifer`` command.  The general syntax for invoking ``pestifer`` is

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
   subs/run
   subs/run-example
   subs/new-system
   subs/fetch-example
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

Pestifer provides a library of useful TcL/VMD procs and VMD macro definitions.  Include the following ``proc`` definition somewhere in your ``.vmdrc`` file.  Then the command ``pestifer_init`` will be available any time you start VMD.  You can issue ``pestifer_init`` from any Tcl prompt in a running VMD session to give it access to pestifer's TcL library.

.. code-block:: tcl

      proc pestifer_init { } {
         set status 0
         if {[catch {exec which pestifer} results options]} {
            set details [dict get $options -errorcode]
            if {[lindex $details 0] eq "CHILDSTATUS"} {
               set status [lindex $details 2]
            } else {
               return -options $options -level 0 $results
            }
         }
         if { $status == 0 } {
            set pestifer_tcl_root [exec pestifer --no-banner wheretcl --root]
            vmdcon -info "Source ${pestifer_tcl_root}/vmdrc.tcl"
            source ${pestifer_tcl_root}/vmdrc.tcl
         } else {
            vmdcon -info "Pestifer is not available in your current environment."
         }
      }

This ``proc`` uses the :ref:`subs wheretcl` subcommand to find the location of pestifer's TcL library, and it sources the ``vmdrc.tcl`` file from that location.  The ``vmdrc.tcl`` file contains the definitions of pestifer's TcL procs, packages, and macros.

If you are using VMD to analyze a system generated using pestifer, it is a good idea to initialize VMD's access to pestifer's TcL library.  This is because pestifer extends the definitions of the ``glycan`` and ``lipid`` macros, among other reasons.  In addition to defining the ``pestifer_init`` proc, you should also add the following line to your ``~/.vmdrc`` file (after the ``pestifer_init`` proc definition):

.. code-block:: tcl

   pestifer_init

This will ensure that the pestifer library is available in your VMD session.  You can then use the ``pestifer`` command to access pestifer's functionality from within VMD scripts.  Note that if you launch VMD from a shell in which the ``pestifer`` command is not available, you will not be able to use the ``pestifer_init`` proc, and you will not have access to pestifer's TcL library.