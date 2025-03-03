.. _example 1:

Example 1: BPTI
---------------

The Input Configuration
=======================

The ``psfgen`` `user manual <https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf>`_ is a great resource for learning how to use ``psfgen``.  A simple example in that manual is a solvation of BPTI starting from its PDB coordinates (PDB ID 6pti).  ``pestifer`` can reproduce this solvation via the input configuration shown below:

.. literalinclude:: ../../../pestifer/resources/examples/01-bpti.yaml
    :language: yaml

You can check the :ref:`config_ref` for a complete reference to Pestifer config files.

This build can be performed (preferably in a clean directory) using this command:

.. code-block:: console

   $ pestifer run-example 1

The first thing ``pestifer`` does with ``run-example`` is to copy the YAML config file for that example into the local directory.  In this case, the file copied is named ``bpti.yaml``, and contains what you see above.  Or, alternatively, pasting that content into a local file ``myconfig.yaml``:

.. code-block:: console

   $ pestifer run myconfig.yaml

You could also use ``fetch-example`` to get the config file and then run it:

.. code-block:: console

  $ ls
  $ pestifer fetch-example 1
  $ ls
  bpti.yaml
  $ pestifer run bpti

(If there is no extension on the argument of run, pestifer assumes one of ``.yaml``, ``.yml``, or ``.ym``.)

``bpti.yaml`` is a YAML-format text file, and the keywords (of course) have particular meanings.  This is also an example of a "minimal" configuration file; ``pestifer`` has many more controls that can be set in a configuration file than are shown here.  Here, this configuration file contains two topmost directives: ``title`` and :ref:`config_ref tasks`.  The value of ``title`` is the string ``BPTI`` and the value of ``tasks`` is a *list*.  Each element in the list of tasks is itself a directive describing a task, and ``pestifer`` in general executes tasks in the order they appear in the ``tasks`` list.

Digression: Interactive Help 
============================

``pestifer`` uses the general-purpose package ``ycleptic`` (`pypi <https://pypi.org/project/ycleptic/>`_) to manage its input configurations.  A package developer using ``ycleptic`` specifies a "pattern" file describing the configuration file syntax they would like their package to have.  ``ycleptic`` provides two useful features:

1. Automatic generaton of a hierarchical arrangement of RST files for documentation of all configuration parameters; in these pages, this is rooted at :ref:`config_ref`.
2. Automatic acquisition of a command-line interactive help feature that allows package users to explore the configuration file format specified by the package developers.  

Let's use this second feature to explore the ``psfgen`` task.  (You can visit the :ref:`config_ref tasks psfgen` page to view the same info in the online documentation.) 

.. code-block:: console

  $ pestifer --no-banner config-help tasks
  Help on user-provided configuration file format

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
  pestifer-help:  psfgen

  psfgen:
      Parameters controlling a psfgen run on an input molecule

  base|tasks->psfgen
      source ->
      mods ->
      cleanup
      .. up
      ! quit
  pestifer-help: source

  source:
      Specifies the source of the initial coordinate file

  base|tasks->psfgen->source
      prebuilt ->
      id
      alphafold
      biological_assembly
      transform_reserves
      remap_chainIDs
      reserialize
      model
      file_format
      cif_residue_map_file
      exclude ->
      sequence ->
      .. up
      ! quit
  pestifer-help:  id

  id:
      The 4-character PDB ID of the source or the basename of a local
        coordinate file (PDB or mmCIF format); pestifer will download
        from the RCSB if a file is not found

This tells us that, in addition to ``id``, we have the ability to set several other control parameters.  Continuing in this interactive help session:

.. code-block:: console

  pestifer-help: biological_assembly

  biological_assembly:
      integer index of the biological assembly to construct; default is 0,
        signifying that the asymmetric unit is to be used
      default: 0

  All subdirectives at the same level as 'biological_assembly':

  base|tasks->psfgen->source
      prebuilt ->
      id
      alphafold
      biological_assembly
      transform_reserves
      remap_chainIDs
      reserialize
      model
      file_format
      cif_residue_map_file
      exclude ->
      sequence ->
      .. up
      ! quit
  pestifer-help: 

And so on.  Let's return to the example.  Immediately after the ``psfgen`` task we declare an ``md`` task, and the subdirective ``ensemble`` is set to ``minimize``.  There are no other subdirectives explicitly listed.  This task will use ``namd3`` to run an energy minimization.  As we did for the ``source`` subdirective of the ``psfgen`` task, let's have a look at the possible subdirectives for an ``md`` task.  We can do this by going "up" twice (``source`` to ``psfgen`` to ``tasks``) and then down into the ``md`` task:

.. code-block:: console

  pestifer-help: ..

  base|tasks->psfgen
      source ->
      mods ->
      cleanup
      .. up
      ! quit
  pestifer-help: ..

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
  pestifer-help: md

  md:
      Parameters controlling a NAMD run

  base|tasks->md
      vacuum
      ensemble
      minimize
      nsteps
      dcdfreq
      xstfreq
      temperature
      pressure
      other_parameters
      constraints ->
      .. up
      ! quit
  pestifer-help:

The Input Configuration (Continued)
===================================

So let's return again to the example.  After this ``md`` task is the ``solvate`` task.  Notice that it has no subdirectives; only default values are used for any subdirectives. Then comes another minimization via an ``md`` task, then an NVT equilibration, and then a series of progressively longer NPT equilibrations in yet more ``md`` tasks.  These "chained-together" NPT runs avoid the common issue that, after solvation, the density of the initial water box is a bit too low, so under pressure control the volume shrinks.  It can shrink so quickly that NAMD's internal data structures for distributing the computational load among processing units becomes invalid, which causes NAMD to die.  The easiest way to reset those internal data structures is just to restart NAMD from the result of the previous run.

The ``mdplot`` task generates a plot of system density (in g/cc) vs time step for the series of MD simulations that occur after solvation.  This is a quick way to check that enough NPT equilibration has been performed.  For this example, the plot looks like this:

.. figure:: pics/solvated-density.png

    Density vs. timestep for the BPTI system post-solvation.

Finally, we see a ``terminate`` task, whose main role is to generate some informative output and to provide a set of NAMD input files (PSF, PDB, xsc, coor, and vel) that all have a common base file name.  The ``package`` subdirective creates a tarball of all required input files to execute a NAMD run, ready for transfer to the HPC resource of your choice.

This run generates a lot of other files.  One such file, ``bpti-complete.yaml`` is the fully explicit configuration file implied by the given configuration file and any default values.  It can be instructive to peruse this file to see the totality of what you can specify for ``pestifer``; it is possible to have very close control over the ``psfgen`` script generation by, for example, adding ``pdbalias`` directives.

The outputs of this build are the PSF/PDB/COOR/VEL/XSC files needed to (re)start namd3; by default, these are ``my_6pti.pdb``, etc.

.. code-block:: console

   $ ls my_6pti*
   my_6pti.coor  my_6pti.pdb  my_6pti.psf  my_6pti.vel  my_6pti.xsc

You should note the presence of CHARMM force-field files in the current directory.  These are generated by ``pestifer`` during the build, and are essentially copies of the parent files with certain lines commented out to permit use by VMD and NAMD.  The parent files are not altered.

.. code-block:: console

  $ tar ztf prod_6pti.tgz
  prod_6pti.namd
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
  02-00-solvate.psf
  08-00-md-NPT.pdb
  08-00-md-NPT.coor
  08-00-md-NPT.xsc
  08-00-md-NPT.vel

``prod_6pti.namd`` is the namd3 configuration file, and it created with some default values.  Carefully consider its contents before you run; you will need to edit it!