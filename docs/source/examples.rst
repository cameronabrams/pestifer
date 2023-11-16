Examples
========

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

This file is YAML format; you can think of it as a python ``dict`` with nesting.  ``pestifer`` uses the general-purpose package ``ycleptic`` (`pypi <https://pypi.org/project/pestifer/>`_) to manage its input configurations.  Under ``ycleptic``, the user provides a YAML format file that contains a set of "directives", where a directive is a ``dict`` with a single key and a value of any type, including directives. Here, there are two topmost directives: ``title`` and ``tasks``.  The value of ``title`` is ``BPTI`` and the value of ``tasks`` is a *list*.  Each element in the list of tasks is itself a directive describing a task, and ``pestifer`` executes tasks in the order they appear in the ``tasks`` list.

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

This tells us that, in addition to ``id``, we have the ability to set seven other keys.  Again, using ``pestifer config-help`` we can learn about these:

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
So let's return again to the example.  After this ``md`` task is the ``solvate`` task.  Notice that it has no subdirectives; only default values are used for any subdirectives. (Currently (v. 1.2.0) the only subdirective is ``pad``.) Then comes another minimization via an ``md`` task, then two sequential NPT equilibrations in yet two more ``md`` tasks.  These "chained-together" NPT runs avoid the common issue that, after solvation, the density of the initial water box is a bit too low, so under pressure control the volume shrinks.  It can shrink so quickly that NAMD's internal data structures for distributing the computational load among processing units becomes invalid, which causes NAMD to die.  The easiest way to reset those internal data structures is just to restart NAMD from the result of the previous run.

Finally, even though it is not explicitly declared here, ``pestifer`` always ends a list of tasks with a special ``terminate`` task, whose main role is to generate some informative output and to provide a set of NAMD input files (PSF, PDB, xsc, coor, and vel) that all have a common base file name.  The default base file name is ``my_system``.

This run will generate a lot of files.  One such file, ``bpti-complete.yaml`` is the fully explicit configuration file implied by the given configuration file and the unstated default values.  It can be instructive to peruse this file to see the totality of what you can specify for ``pestifer``; it is possible to have very close control over the ``psfgen`` script generation by, for example, adding ``pdbalias`` directives.

The outputs of this build are the PSF/PDB/COOR/VEL/XSC files needed to (re)start namd2; by default, these are ``my_system.pdb`` etc.

.. code-block:: console

   $ ls my_system*
   my_system.coor  my_system.pdb  my_system.psf  my_system.vel  my_system.xsc

Where, you may wonder, are the CHARMM parameter files?  ``pestifer`` includes the full CHARMM force field download from the MacKerrel lab, so they are somewhere under your python package installation tree.  You can see their full pathnames in any NAMD config file ``pestifer`` generates along the way.  However, to make things easier, we include a ``package`` subdirective in the ``terminate`` task that will copy the necessary CHARMM parameter files to the local directory for you.  We illustrate this in Example 11.

Example 11
----------
This is the same as Example 1, except we delete the phosphate ion, and we introduce a ``terminate`` task with a ``package`` directive.  ``package`` will make a tarball of all necessary input files for a production ``namd2`` run, including the PSF/PDB/COOR/VEL/XSC fileset and all parameter files.

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

Note the ``exclude`` subdirective under ``source``.  You remember how you can learn about it?  Using ``config-help``: 

.. code-block:: console

  $ pestifer config-help tasks psfgen source exclude --no-banner
  Help on user-provided configuration file format
  tasks->
  psfgen->
  source->
  exclude:
      Specifies any residues or atoms present in the PDB source to exclude
        from the system
      type: dict
      Help available for chains, resnames

Each of ``chains`` and ``resnames`` are lists, and in the configuration file above, we have a single-element list for ``resnames`` that indicates the resname ``PO4``, which is how the phosphate ion is labelled in the original PDB file.

After the build, the name of the tarball generated is ``prod_6pti.tgz`` and its contents are:

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

``prod_6pti.namd`` is the NAMD2 configuration file, set up with some default values.  Carefully consider its contents before you run!

