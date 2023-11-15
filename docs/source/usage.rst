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

Example 1
---------

An example YAML configuration file to build the 6PTI system described in the 
in the `psfgen user manual <https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf>`_ is below:

.. code-block:: yaml

  title: BPTI
  tasks:
    - psfgen:
        source:
          id: 6pti
    - md:
        ensemble: minimize
    - solvate:
    - md:
        ensemble: minimize
    - md:
        ensemble: NPT
        nsteps: 1000
    - md:
        ensemble: NPT
        nsteps: 1000

This build can be performed (preferably in a clean directory) using this command:

.. code-block:: console

   $ pestifer run-example 1

Or, alternatively, pasting that content into a local file ``myconfig.yaml``:

.. code-block:: console

   $ pestifer run myconfig.yaml

The first thing ``pestifer`` does with ``run-example`` is to copy the YAML config file for that example into the local directory.  In this case, the file copied is named ``bpti.yaml``, and contains what you see above.

This file is YAML format; you can think of it as a python ``dict`` with nesting.  ``pestifer`` uses the general-purpose package ``ycleptic`` (`pipy <https://pypi.org/project/pestifer/>`_) to manage its input configurations.  Under ``ycleptic``, the user provides a YAML format file that contains a set of "directives", where a directive is a ``dict`` with a single key and a value of any type, including directives. Here, there are two topmost directives: ``title`` and ``tasks``.  The value of ``title`` is ``BPTI`` and the value of ``tasks`` is a *list*.  Each element in the list of tasks is itself a directive describing a task, and ``pestifer`` executes tasks in the order they appear in the ``tasks`` list.

For the ``psfgen`` task, we see the directive ``source``.  Its value appears to be yet another subdirective, ``id``, but the value of source is a ``dict`` with several keys, and we specify *only* ``id``, and the others are set to default values.  We can see these other keys and their default values using ``pestifer config-help``: 

.. code-block:: console

  $ pestifer config-help tasks psfgen source --no-banner
  Help on user-provided configuration file format
  tasks->
  psfgen->
  source:
      Specifies the source of the initial coordinate file
      type: dict
      Help available for id, biological_assembly, file_format, cif_residue_map_file, psf, altcoords, exclude, sequence

This tells us that, in addition to `id`, we have the ability to set seven other keys.  Again, using `pestifer config-help` we can learn about these:

.. code-block:: console

  $ pestifer config-help tasks psfgen source id --no-banner
  Help on user-provided configuration file format
  tasks->
  psfgen->
  source->
  id:
      The 4-character PDB ID of the source or the basename of a local
        coordinate file (PDB or mmCIF format); pestifer will download
        from the RCSB if a file is not found
      type: str
      A value is required.
  $ pestifer config-help tasks psfgen source biological_assembly --no-banner
  Help on user-provided configuration file format
  tasks->
  psfgen->
  source->
  biological_assembly:
      integer index of the biological assembly to construct; default is 0,
        signifying that the asymmetric unit is to be used
      type: int
      default: 0
  $ pestifer config-help tasks psfgen source file_format --no-banner
  Help on user-provided configuration file format
  tasks->
  psfgen->
  source->
  file_format:
      either PDB or mmCIF; some entries do not have a PDB-format file.  The
        main advantage of PDB is that it uses the author-designations
        for chains by default.  mmCIF is the new "default" format of the
        PDB.
      type: str
      default: PDB
      allowed values: PDB, mmCIF

And so on.  Let's return to the example.  Immediately after the ``psfgen`` task we declare an ``md`` task, and the subdirective ``ensemble`` is set to ``minimize``.  There are no other subdirectives explicitly listed.  This task will use ``namd2`` to run an energy minimization.  As we did for the ``source`` subdirective of the ``psfgen`` task, let's have a look at the possible subdirectives for an ``md`` task:

.. code-block:: console

  $ pestifer config-help tasks md --no-banner
  Help on user-provided configuration file format
  tasks->
  md:
      Parameters controlling a NAMD run
      type: dict
      Help available for ensemble, minimize, nsteps, dcdfreq, xstfreq, temperature, pressure, other_parameters, constraints

By now, you know how to use ``config-help`` to figure out what these subdirectives mean. 
So let's return again to the example.  After this ``md`` task is the ``solvate`` task.  Notice that it has _no_ subdirectives; only default values are used for any subdirectives. (Currently (v. 1.2.0) the only subdirective is ``pad``.) Then comes another minimization via an ``md`` task, then two sequential NPT equilibrations in yet two more ``md`` tasks.  These "chained-together" NPT runs avoid the common issue that, after solvation, the density of the initial water box is a bit too low, so under pressure control the volume shrinks.  It can shrink so quickly that NAMD's internal data structures for distributing the computational load among processing units becomes invalid, which causes NAMD to die.  The easiest way to reset those internal data structures is just to restart NAMD from the result of the previous run.

Finally, even though it is not explicitly declared here, ``pestifer`` always ends a list of tasks with a special ``terminate`` task, whose main role is to generate some informative output and to provide a set of NAMD input files (PSF, PDB, xsc, coor, and vel) that all have a common base file name.  The default base file name is ``my_system``.

This run will generate a lot of files.  One such file, ``bpti-complete.yaml`` is the fully explicit configuration file implied by the given configuration file and the unstated default values.  It can be instructive to peruse this file to see the totality of what you can specify for ``pestifer``; it is possible to have very close control over the ``psfgen`` script generation by, for example, adding ``pdbalias`` directives.

The outputs of this build are the PSF/PDB/COOR/VEL/XSC files needed to (re)start namd2; by default, these are ``my_system.pdb`` etc.

.. code-block:: console

   $ ls my_system*
   my_system.coor  my_system.pdb  my_system.psf  my_system.vel  my_system.xsc

Example 2
---------
This is the same as Example 1, except we specify that the phosphate ion present in the PDB input file be excluded.

.. code-block:: yaml

  title: BPTI with phosphate ion excluded
  tasks:
    - psfgen:
        source:
          id: 6pti
          exclude:
            resnames:
              - PO4
    - md:
        ensemble: minimize
    - solvate:
    - md:
        ensemble: minimize
    - md:
        ensemble: NPT
        nsteps: 1000
    - md:
        ensemble: NPT
        nsteps: 1000

This exclusion is specified via the ``exclude`` directive under ``source``, and we further specify that it is all residues with a particular name that we are excluding (``PO4``).  We must refer to any resnames to exclude using exactly the same string that refers to them in the PDB input itself.

Example 11
----------
This is the same as Example 2, except we introduce a ``terminate`` task with a ``package`` directive.  ``package`` will make a tarball of all necessary input files for a production ``namd2`` run, including the PSF/PDB/COOR/VEL/XSC fileset and all parameter files.

.. code-block:: yaml

  title: BPTI, packaging all inputs for NAMD deployment
  tasks:
    - psfgen:
        source:
          id: 6pti
          exclude:
            resnames:
              - PO4
    - md:
        ensemble: minimize
    - solvate:
    - md:
        ensemble: minimize
    - md:
        ensemble: NVT
    - md:
        ensemble: NPT
        nsteps: 500
    - md:
        ensemble: NPT
        nsteps: 1500
    - terminate:
        basename: my_6pti
        package:
          ensemble: NPT
          basename: prod_6pti

The name of the tarball generated is ``prod_6pti.tgz`` and its contents are:

.. code-block:: console

   $ tar ztf prod_6pti.tgz
   par_all36m_prot.prm
   par_all36_carb.prm
   par_all36_lipid.prm
   par_all36_carb.prm
   par_all36_na.prm
   par_all36_cgenff.prm
   toppar_all36_carb_glycopeptide.str
   toppar_all36_prot_modify_res.str
   toppar_water_ions.str
   toppar_all36_moreions.str
   01-00-solvate.psf
   02-00-relax-relax.pdb
   02-00-relax-relax.coor
   02-00-relax-relax.xsc
   02-00-relax-relax.vel
   prod_6pti.namd
   $

Carefully consider the contents of the ``namd2`` config file before you run!

