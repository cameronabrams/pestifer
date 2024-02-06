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

There are several subcommands:

``run``: the main subcommand that uses a user's input config file to build a system.

.. code-block:: console

   $ pestifer run <config.yaml>

Here ``config.yaml`` is the name of the configuration file that describes the build.  The
best way to learn about the configuration file is to run the examples provided, and
then you can try to make your own.

``run-example``: there are 12 example systems; to run number four, for example:

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

In VMD scripts
--------------

Pestifer has a pretty handy library of TcL procedures.  If you want to peruse the source, pestifer will tell you where to find them:

.. code-block:: console

   $ ls `pestifer wheretcl --proc-dir`
   autools.tcl  checkpierce.tcl  crot.tcl  dcdlog.tcl  declash.tcl  getlinks.tcl  multimer.tcl  numbering.tcl  saverestore.tcl  util.tcl

If you want to use any of the procs defined in those files in your own VMD script, the easiest thing to do is to put this proc definition in your own VMD startup file:

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
         set script_dir [exec pestifer wheretcl --script-dir]
         vmdcon -info "Sourcing ${script_dir}/pestifer-vmd.tcl"
         return ${script_dir}/pestifer-vmd.tcl
      } else {
         vmdcon -info "Pestifer is not available in your current environment."
      }
   }

Then, you can use it in a source command in any VMD script or TcL session you like:

.. code-block:: tcl

   source [pestifer_init]

