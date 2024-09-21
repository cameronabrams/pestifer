Example 1: BPTI
---------------

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
        ensemble: NVT
    - md:
        ensemble: NPT
        nsteps: 200
    - md:
        ensemble: NPT
        nsteps: 400
    - md:
        ensemble: NPT
        nsteps: 800
    - md:
        ensemble: NPT
        nsteps: 1600
    - mdplot:
        savedata: solvated.csv
        traces:
          - density
        units:
          density: g_per_cc
        basename: solvated        
    - terminate:
        basename: my_6pti
        package:
          ensemble: NPT
          basename: prod_6pti

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

This file is YAML format; you can think of it as a python ``dict`` with nesting.  ``pestifer`` uses the general-purpose package ``ycleptic`` (`pypi <https://pypi.org/project/ycleptic/>`_) to manage its input configurations.  Under ``ycleptic``, the user provides a YAML format file that contains a set of "directives", where a directive is a ``dict`` with a single key and a value of any type, including directives. Here, there are two topmost directives: ``title`` and ``tasks``.  The value of ``title`` is the string ``BPTI`` and the value of ``tasks`` is a *list*.  Each element in the list of tasks is itself a directive describing a task, and ``pestifer`` executes tasks in the order they appear in the ``tasks`` list.

For the ``psfgen`` task, we see the directive ``source``.  Its value appears to be yet another subdirective, ``id``, but the value of source is a ``dict`` with several keys, and we specify *only* ``id``, and the others are set to default values.  We can see these other keys and their default values using ``pestifer config-help``: 

.. code-block:: console

  $ pestifer --no-banner config-help tasks psfgen source
  Help on user-provided configuration file format
  tasks->
  psfgen->
  source:
      Specifies the source of the initial coordinate file
      type: dict
      Help available for id, biological_assembly, file_format, cif_residue_map_file, psf, altcoords, exclude, sequence

This tells us that, in addition to ``id``, we have the ability to set seven other keys.  Again, using ``pestifer config-help`` we can learn about these:

.. code-block:: console

  $ pestifer --no-banner config-help tasks psfgen source id 
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
  $ pestifer --no-banner config-help tasks psfgen source biological_assembly
  Help on user-provided configuration file format
  tasks->
  psfgen->
  source->
  biological_assembly:
      integer index of the biological assembly to construct; default is 0,
        signifying that the asymmetric unit is to be used
      type: int
      default: 0
  $ pestifer --no-banner config-help tasks psfgen source file_format
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

  $ pestifer --no-banner config-help tasks md
  Help on user-provided configuration file format
  tasks->
  md:
      Parameters controlling a NAMD run
      type: dict
      Help available for ensemble, minimize, nsteps, dcdfreq, xstfreq, temperature, pressure, other_parameters, constraints

By now, you know how to use ``config-help`` to figure out what these subdirectives mean. 
So let's return again to the example.  After this ``md`` task is the ``solvate`` task.  Notice that it has no subdirectives; only default values are used for any subdirectives. Then comes another minimization via an ``md`` task, then an NVT equilibration, and then a series of progressively longer NPT equilibrations in yet more ``md`` tasks.  These "chained-together" NPT runs avoid the common issue that, after solvation, the density of the initial water box is a bit too low, so under pressure control the volume shrinks.  It can shrink so quickly that NAMD's internal data structures for distributing the computational load among processing units becomes invalid, which causes NAMD to die.  The easiest way to reset those internal data structures is just to restart NAMD from the result of the previous run.

The ``mdplot`` task generates a plot of system density (in g/cc) vs time step for the series of MD simulations that occur after solvation.  This is a quick way to check that enough NPT equilibration has been performed.  For this example, the plot looks like this:

.. figure:: pics/solvated-density.png

    Density vs. timestep for the BPTI system post-solvation.

Finally, we see a ``terminate`` task, whose main role is to generate some informative output and to provide a set of NAMD input files (PSF, PDB, xsc, coor, and vel) that all have a common base file name.  The ``package`` subdirective creates a tarball of all required input files to execute a NAMD run, ready for transfer to the HPC resource of your choice.

This run generates a lot of other files.  One such file, ``bpti-complete.yaml`` is the fully explicit configuration file implied by the given configuration file and any default values.  It can be instructive to peruse this file to see the totality of what you can specify for ``pestifer``; it is possible to have very close control over the ``psfgen`` script generation by, for example, adding ``pdbalias`` directives.

The outputs of this build are the PSF/PDB/COOR/VEL/XSC files needed to (re)start namd2; by default, these are ``my_6pti.pdb``, etc.

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

``prod_6pti.namd`` is the NAMD2 configuration file, and it created with some default values.  Carefully consider its contents before you run; you will need to edit it!