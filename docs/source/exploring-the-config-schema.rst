.. _exploring_config_schema:

Exploring the configuration schema
==================================


``pestifer`` uses the general-purpose package ``ycleptic`` (`pypi <https://pypi.org/project/ycleptic/>`_) to manage its input configurations.  A package developer using ``ycleptic`` specifies a "pattern" file describing the configuration file syntax they would like their package to have.  ``ycleptic`` provides two useful features:

1. Automatic generaton of a hierarchical arrangement of RST files for documentation of all configuration parameters; in these pages, this is rooted at :ref:`config_ref`.
2. Automatic acquisition of a command-line interactive help feature that allows package users to explore the configuration file format specified by the package developers.  

Let's use this second feature to explore the ``fetch`` task.  (You can visit the :ref:`config_ref tasks` page to view the same info in the online documentation.) 

.. code-block:: bash

  $ pestifer --no-banner config-help tasks

    tasks:
        Specifies the tasks to be performed serially in a pestifer run

    base|tasks
        fetch ->
        continuation ->
        merge ->
        psfgen ->
        ligate ->
        pdb2pqr ->
        mdplot ->
        cleave ->
        solvate ->
        desolvate ->
        ring_check ->
        make_membrane_system ->
        md ->
        density_equilibrate ->
        membrane_equilibrate ->
        manipulate ->
        terminate ->
        validate ->
        .. up
        ! quit
    pestifer-help: fetch

    fetch:
        Fetch task; its only job is to fetch any external data file (e.g.,
          PDB).

    base|tasks->fetch
        source
        sourceID
        source_format
        .. up
        ! quit
    pestifer-help: source

    source:
        Source for the initial coordinate file; one of ``rcsb`` (for the RCSB
          PDB), ``alphafold`` (for the AlphaFold DB), ``opm`` (for the OPM
          database), or ``local`` (for a local file)
        default: rcsb
        allowed values: rcsb, alphafold, opm, local

    All subattributes at the same level as 'source':

    base|tasks->fetch
        source
        sourceID
        source_format
        .. up
        ! quit
    pestifer-help: sourceID

    sourceID:
        ID of the source file; if source is ``local``, a file ``sourceID.pdb``
          or ``sourceID.cif`` must exist in the working directory

    All subattributes at the same level as 'sourceID':

    base|tasks->fetch
        source
        sourceID
        source_format
        .. up
        ! quit
    pestifer-help: source_format

    source_format:
        Format of the source file; this should be ``pdb`` or ``cif``
        default: pdb
        allowed values: pdb, cif

    All subattributes at the same level as 'source_format':

    base|tasks->fetch
        source
        sourceID
        source_format
        .. up
        ! quit
    pestifer-help: !
  $

In the config file for this example, we specify on the the ``sourceID`` as 6pti; the other source attributes take their default values.  This causes ``pestifer`` to fetch the file ``6pti.pdb`` from the RCSB PDB (if ``6pti.pdb`` does not already exist in the current working directory).

We can return to ``config-help`` to explore the ``psfgen`` task, which is the next task in the list.  We can do this by:

.. code-block:: bash

  $ pestifer config-help tasks psfgen

    psfgen:
        Parameters controlling a specific psfgen run on an input molecule

    base|tasks->psfgen
        source ->
        mods ->
        .. up
        ! quit
    pestifer-help: source

    source:
        Specifies the processing and interpretation of the initial source
          coordinate file

    base|tasks->psfgen->source
        biological_assembly
        transform_reserves
        remap_chainIDs
        reserialize
        model
        cif_residue_map_file
        include
        exclude
        sequence ->
        .. up
        ! quit
    pestifer-help: biological_assembly

    biological_assembly:
        integer index of the biological assembly to construct; default is 0,
          signifying that the asymmetric unit is to be used
        default: 0

    All subattributes at the same level as 'biological_assembly':

    base|tasks->psfgen->source
        biological_assembly
        transform_reserves
        remap_chainIDs
        reserialize
        model
        cif_residue_map_file
        include
        exclude
        sequence ->
        .. up
        ! quit
    pestifer-help:

The same works for any task.  The ``md`` task, for instance, carries a long list of
subdirectives, none of which you have to memorize:

.. code-block:: bash

  $ pestifer config-help tasks md

    md:
        Parameters controlling a NAMD run

    base|tasks->md
        single-core
        cpu-override
        vacuum
        ensemble
        minimize
        nsteps
        dcdfreq
        xstfreq
        temperature
        pressure
        addl_paramfiles
        other_parameters
        constraints ->
        ssrestraints ->
        colvar_specs
        .. up
        ! quit
    pestifer-help:

Every one of these is also in the :ref:`config_ref`, which is generated from the same
schema, so the interactive explorer and the online reference can never disagree.
