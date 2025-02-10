config-help
-----------

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
