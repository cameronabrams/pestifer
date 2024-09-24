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
...........

``run`` 
+++++++

Builds a system according to the input configuration file.

.. code-block:: console

   $ pestifer run <config.yaml>

Here ``config.yaml`` is the name of the configuration file that describes the build.  One way to learn about the configuration file is to run the examples provided, and
then you can try to make your own.

``run-example``
+++++++++++++++

There are 14 example systems; to run number four, for example:

.. code-block:: console
   
   $ pestifer run-example 4

(Best to do that in a clean directory.)  

``fetch-example``
+++++++++++++++++

Like ``run-example``, except it only copies the YAML input file needed to run the example to the CWD.  It can be edited, if desired, and then run using the ``run`` subcommand (see above).  If you request help on ``fetch-example``, you will see a list of examples:

.. code-block:: console

   $ pestifer fetch-example --help
   usage: pestifer fetch-example [-h] [--config-updates CONFIG_UPDATES [CONFIG_UPDATES ...]] number

   Fetch YAML config for one of the examples:
   1: BPTI
   2: BPTI with phosphate ion excluded
   3: BPTI, no phosphate, some random mutations plus deletion of one disulfide
   4: BPTI, no phosphate, introducing a disulfide via mutations
   5: HIV-1 Env Trimer 4zmj
   6: HIV-1 Env Trimer 4tvp, Fabs and sulfate ions removed
   7: HIV-1 Env Trimer 8fad, drug molecule removed
   8: HIV-1 Env Trimer 8fae, drug molecule removed
   9: HIV-1 Env Trimer 7txd, substitutions, no ligands
   10: HIV-1 Env Trimer 5vn3, excluding sCD4 and Fab chains, with Gly3 stubs replacing missing v1/2
   11: hexameric insulin
   12: insulin receptor ectodomain, Fabs removed
   13: BA.2 SARS-CoV-2 Spike 7xix, fully glycosylated using grafts, and cleaved
   14: HIV-1 gp41 MPER-TM trimer 6e8w embedded in a model viral membrane

   positional arguments:
   number                example number

   options:
   -h, --help            show this help message and exit
   --config-updates CONFIG_UPDATES [CONFIG_UPDATES ...]
                           yaml files to update example

``config-help``
+++++++++++++++

Interactive help in constructing a config file. Because it uses ``ycleptic``, ``pestifer`` has a built-in interactive system for help generating YAML-format input configuration files.  

.. code-block:: console

   $ pestifer config-help
   Help on user-provided configuration file format
      charmmff ->
      psfgen ->
      namd2 ->
      title
      paths ->
      tasks ->
      .. up
      ! quit
   pestifer-help: 

Each of these topics is a top-level directive allowed in a config file.  We can dig down on any one of them:

.. code-block:: console

   pestifer-help: tasks

   tasks:
      Specifies the tasks to be performed serially

   base|tasks
      restart ->
      psfgen ->
      ligate ->
      mdplot ->
      cleave ->
      domainswap ->
      solvate ->
      ring_check ->
      bilayer ->
      md ->
      manipulate ->
      terminate ->
      .. up
      ! quit
   pestifer-help:

Continuing to drill down is easy -- just add the next directive to the interactive-help command line:

.. code-block:: console

   pestifer-help: bilayer

   bilayer:
      Parameters controlling packmol to generate a bilayer system

   base|tasks->bilayer
      embed ->
      lipids
      mole_fractions
      solvents
      solvent_mole_fractions
      SAPL
      leaflet_thickness
      scale_excluded_volume
      fuzz_factor
      dims
      length_pad
      solution_gcc
      cation
      anion
      salt_con
      nloop
      nloop_all
      tolerance
      seed
      .. up
      ! quit
   pestifer-help: 

Using Pestifer in your own VMD scripts and sessions
---------------------------------------------------

Pestifer has a pretty handy library of TcL packages.  If you want to peruse the source, pestifer will tell you where to find them:

.. code-block:: console

   $ pestifer wheretcl --pkg-dir

If you want to use any of the procs defined in those packages in your own VMD script, the easiest thing to do is to put this ``proc`` definition in your own VMD startup file:

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
            return ${pestifer_tcl_root}/vmdrc.tcl
         } else {
            vmdcon -info "Pestifer is not available in your current environment."
         }
      }

Then, you can use it in a source command in any VMD script or TcL session you like:

.. code-block:: tcl

   source [pestifer_init]

This of course requires that your VMD session was launched from a shell running a python virtual environment in which ``pestifer`` is installed.
