.. _subs_runtasks_psfgen:

psfgen 
------

A ``psfgen`` task is usually the first task in a ``tasks`` directive.  Its basic functionality is to set up and conduct the first ``psfgen`` run to generate psf and pdb files from a structure file.  It's real job is writing the ``psfgen`` input file and then calling ``psfgen``.

An example ``psfgen`` task that fetches and processes the 6PTI structure might look like this:

.. code-block:: yaml

  tasks:
    - psfgen:
        source:
          id: 6pti  # if 6pti.pdb does not exist in the current directory, it will be fetched from the RCSB
    - ... (more tasks follow)

With this input, pestifer will fetch ``6pti.pdb`` from the RCSB, generate the required ``psfgen`` script, and execute ``psfgen`` from within a text-only VMD session.  The PSF and PDB files that result are passed forward to the next task in the ``tasks`` directive.  Any string is allowed in an ``id`` directive, but if ``<id-value>.pdb`` does not exist in the current directory, pestifer will attempt to fetch it from the RCSB PDB.  If you want to use a local PDB file, you can simply place it in the current directory and use the ``id`` directive to specify its name (without the ``.pdb`` extension).

The only two subdirectives under ``psfgen`` are ``source`` and ``mods``. 

source
++++++

Under ``source``, there are three mutually exclusive directives used to specify the molecular input:

1. ``id``: This is a 4-character PDB ID.
2. ``alphafold``: This is an AlphaFold database ID.
3. ``prebuilt``: This allows you to specify obligatory PSF, PDB, and optional xsc file for a system that was already built.

Secondary directives under ``source`` are described in its :ref:`Configuration File Reference pages <config_ref tasks psfgen source>`.

If ``psfgen`` is the first task, the name of the ``psfgen`` input file it will write is ``00-00-00_psfgen-build.tcl`` and the ``psfgen`` session log will be written to ``00-00-00_psfgen-build.log``. (File naming conventions are briefly explained :ref:`here: <file name conventions>`.  There is no ``pestifer`` use case in which is top-level ``psfgen`` task is *not* the first task in a ``tasks`` directive.)  

A common way of iterating with pestifer when a build is not successful is to carefully read the log file, and then edit the ``00-00-00_psfgen-build.tcl`` file to attempt fix the problem. For example, you might need to add a missing patch or modify a residue name. This requires a good working understanding of how to use ``psfgen``, of course.  You can then re-run the ``psfgen`` task by running the following command in a text-only VMD session:

.. code-block:: bash

   $ vmd -dispdev text
   VMD> pestifer_init
   VMD> source 00-00-00_psfgen-build.tcl 

Once you have a successful standalone ``psfgen`` run, you can then use the ``pestifer`` command to complete the build by reading in the successfully (and manually) created PSF and PDB files in a ``prebuilt`` directive, rather then ``id`` or ``alphafold``.  Or, you may figure out how to modify your original config file to fix the problem and then re-run the ``pestifer`` command with the new config file.  This is a common workflow when using pestifer to build systems, and it is preferred if you want to use pestifer in a high-throughput setting.

``pestifer_init`` is a TcL proc that you should put in your ``~/.vmdrc`` file.  It sets up the environment for pestifer within the VMD session. (See :ref:`use in vmd scripts` for details on how to set up your ``~/.vmdrc`` file to make sure VMD works with pestifer.)

mods 
++++

Under ``mods``, you can specify any number of modifications that can be handled in a single invocation of ``psfgen``.  These include point mutations, deletions, substitutions, additional disulfide bonds or other covalent bonds, grafting of complete glycans onto glycan stems, among others.  

For example, a ``psfgen`` task that builds a 6PTI system with two specific point mutations might look like this:

.. code-block:: yaml

  tasks:
    - psfgen:
        source:
          id: 6pti
        mods:
          mutations:
            - A:T11A # threonine to alanine at position 11, 
            - A:PRO,13,ALA # alternate syntax; proline to alanine at position 13

Below are all the ``mods`` that can be invoked in a ``mods`` directive.  Each of these is described in detail in its :ref:`Configuration File Reference pages <config_ref tasks psfgen mods>`.

.. toctree::
   :maxdepth: 1

   mods/patches
   mods/mutations
   mods/deletions
   mods/insertions
   mods/crotations
   mods/transrot
   mods/ssbonds
   mods/links
   mods/substitutions
   mods/grafts
