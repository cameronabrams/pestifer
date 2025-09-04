.. _subs_runtasks_psfgen:

psfgen 
------

A ``psfgen`` task's  basic functionality is to set up and conduct the first ``psfgen`` run to generate psf and pdb files from a structure file.  It's real job is writing the ``psfgen`` input script and then executing VMD with that script as input.

A traditional ``psfgen`` workflow using VMD begins with the preprocessing the downloaded structure file (if necessary), and then writing a running a ``psfgen`` script.  Pestifer's ``psfgen`` task handles all preprocessing, and additionally can perform several types of modifications to the source structure on its way to creating simulation-ready PSF and PDB files.

An example Pestifer input that specifies fetching and processing the 6PTI structure might look like this:

.. code-block:: yaml

  tasks:
    - fetch:
        sourceID: 6pti
    - psfgen:
    - ... (more tasks follow)

With this input, pestifer will fetch ``6pti.pdb`` from the RCSB, generate the required ``psfgen`` script, and execute VMD with that script as input.  The PSF and PDB files that result are passed forward to the next task in the ``tasks`` directive.

A common way of iterating with pestifer when a build is not successful is to carefully read the log file, and then edit the ``00-00-00_psfgen-build.tcl`` file to attempt fix the problem. For example, you might need to add a missing patch or modify a residue name. This requires a good working understanding of how to use ``psfgen``, of course.  You can then re-run the ``psfgen`` task by running the following command in a text-only VMD session:

.. code-block:: bash

   $ vmd -dispdev text
   VMD> pestifer_init
   VMD> source 00-00-00_psfgen-build.tcl 

Once you have a successful standalone ``psfgen`` run, you can then use the ``pestifer`` command to complete the build by reading in the successfully (and manually) created PSF and PDB files in a ``continuation`` task.  Or, you may figure out how to modify your original config file to fix the problem and then re-run the ``pestifer`` command with the new config file.  This is a common workflow when using pestifer to build systems, and it is preferred if you want to use pestifer in a high-throughput setting.

``pestifer_init`` is a TcL proc that you should put in your ``~/.vmdrc`` file.  It sets up the environment for pestifer within the VMD session. (See :ref:`use in vmd scripts` for details on how to set up your ``~/.vmdrc`` file to make sure VMD works with pestifer.)

The only two subdirectives under ``psfgen`` are ``source`` and ``mods``. 

.. _subs_runtasks_psfgen_source:

source
++++++

Under ``source``, there are several useful directives:

1. ``biological_assembly``: The integer ID of the biological assembly represented in the input structure file which you would like to use to prepare your system.  This is useful for structures that have multiple biological assemblies.  If not specified, the asymmetric unit will be used.
2. ``transform_reserves``: A dictionary mapping chainID's in the asymmetric unit to their counterparts in other protomers of the biological assembly (if the assembly replicates the asymmetric unit).
3. ``remap_chainIDs``: A dictionary mapping chainID's in the input structure file to their counterparts in the output structure file.
4. ``reserialize``: A boolean. If true, resets the serial numbers of all atoms based on their order in the input structure file -- use with caution!!
5. ``model``: The integer ID of the model to use from a multi-model structure file.  If unspecified, the **last** model in the file is used.
6. ``cif_residue_map_file``: If the structure file's format is CIF (rather than PDB), it often has its own labels for chainIDs and residue sequence numbers that differ from what the authors use in their publications.  This specifies the name of an output file that contains the information mapping the CIF-format resids to the author-format.
7. ``include``: A list of strings representing **pythonic** logic for including atoms in the structure file in the final build.
8. ``exclude``: A list of strings representing **pythonic** logic for excluding atoms from the structure file in the final build.
9. ``sequence``: A dictionary with special keyword-value pairs for handling certain sequence modifications to the base molecule. These are:

     a. ``include_terminal_loops``: A boolean.  If true, any unresolved residues in the N and C termini of **all** chains in the input strucutre file will be grown in. Default is **false**.
     b. ``build_zero_occupancy_C_termini``: A list of chainIDs for which you would specifically like to build in any unresolved C-terminal residues.  Default is **empty**.
     c. ``build_zero_occupancy_N_termini``: A list of chainIDs for which you would specifically like to build in any unresolved N-terminal residues.  Default is **empty**.
     d. ``fix_engineered_mutations``: A boolean.  If true, any engineered mutations specified in the input structure file will be reverted back to wild-type in the output structure file.  Default is **false**.
     e. ``fix_conflicts``: A boolean.  If true, any sequence conflicts listed in the input structure file are resolved in favor of the database values.  Default is **false**.
     f. ``loops``: A dictionary with parameters governing how missing loops are grown in.
     g. ``glycans``: A dictionary with parameters governing how glycans are modeled.

source.sequence.loops
#####################

Under the ``loops`` directive under ``source.sequence``, one can specify three main directives.

1. ``sac_res_name``: The is the 3-letter name of the residue used as a temporary C-terminus on the growing loop.  This is necessary to separate the actual least residue of the loop from an active C-terminus, making it very easy to delete the C-terminus prior to linking the C-terminus to the next residue to close the gap.  Defaul is **GLY**.
2. ``min_loop_length``: The minimum number of residues in a loop that qualifies it for the steering simulations to close the gap.  Default is **4**.
3. ``declash``: A dictionary of parameters governing a declashing procedure for grown in loop prior to steering.

    a. ``maxcycles``: maximum number of declash cycles; a cycle is a random torsion angle displacment per residue of the loop.  Default is **0**, which turns off declashing.
    b. ``include_C_termini``: A boolean.  If true, any C-terminal loops are subject to declashing (provided that ``maxcycles`` is greater than 0.)  Default is **true**.
    c. ``clash_dist``: The minimum distance in Angstrom between any two atoms in the loop during the declashing procedure.  Default is **1.5**.


source.sequence.glycans
#######################

Under the ``glycans`` directive under ``source.sequence``, one can specify a single ``declash`` directive, which like the loop declashing, has ``maxcycles`` (default **0**) and ``clashdist`` (default **1.5**).


mods 
++++

Under ``mods``, you can specify any number of modifications that can be handled in a single invocation of ``psfgen``.  These include point mutations, deletions, substitutions, additional disulfide bonds or other covalent bonds, grafting of complete glycans onto glycan stems, among others.  

For example, a ``psfgen`` task that builds a 6PTI system with two specific point mutations might look like this:

.. code-block:: yaml

  tasks:
    - fetch:
        sourceID: 6pti
     - psfgen:
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
