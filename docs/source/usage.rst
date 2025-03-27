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
   subs/cleanup

Use in VMD scripts
------------------

Pestifer provides a library of useful TcL/VMD procs and VMD macro definitions.  If you define the following proc in your `.vmdrc` file, you can access the pestifer library from any VMD script or TcL session.  This proc will automatically source the pestifer library when you start VMD.

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

Then, you can use it in a source command in any VMD script or TcL session you like: 

.. code-block:: tcl

   pestifer_init

If you are using VMD to analyze a system generated using pestifer, it is a good idea to initialize VMD's access to pestifer's TcL library.