config-help
-----------

Because it uses ``ycleptic``, ``pestifer`` has a built-in, interactive, command-line system for help generating YAML-format input configuration files available through the ``config-help`` subcommand.  

.. code-block:: bash

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

This command ends at a prompt (``pestifer-help:``) that allows you to drill down into the help system.  The help system is organized around the topics that are allowed in a config file, which appear in the list above.  Any item with an arrow after it can be drilled down into.  Double-dot (``..``) takes you up, and bang (``!``) quits.  For example, in the case above, if you type ``tasks``, you will get a list of the tasks that can be performed in a config file:

.. code-block:: bash

   pestifer-help: tasks

   tasks:
       Specifies the tasks to be performed serially in a pestifer run

   base|tasks
       restart ->
       psfgen ->
       ligate ->
       mdplot ->
       cleave ->
       domainswap ->
       solvate ->
       desolvate ->
       ring_check ->
       make_membrane_system ->
       md ->
       manipulate ->
       terminate ->
       .. up
       ! quit
   pestifer-help:

Continuing to drill down is easy -- just add the next directive to the interactive-help command line:

.. code-block:: bash

   pestifer-help: make_membrane_system

   make_membrane_system:
       Parameters controlling packmol to generate a membrane system with an
         embedded protein

   base|tasks->make_membrane_system
       bilayer ->
       embed ->
       .. up
       ! quit
   pestifer-help: bilayer

   bilayer:
      Parameters controlling bilayer generation

   base|tasks->make_membrane_system->bilayer
      prebuilt ->
      lipids
      lipid_conformers
      mole_fractions
      patch_nlipids ->
      composition ->
      half_mid_zgap
      solvents
      solvent_mole_fractions
      solvent_to_lipid_ratio
      SAPL
      dims
      npatch
      solution_gcc
      cation
      anion
      salt_con
      nloop
      nloop_all
      tolerance
      seed
      relaxation_protocols ->
      .. up
      ! quit