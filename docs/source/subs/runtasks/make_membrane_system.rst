.. _subs_runtasks_make_membrane_system:

make_membrane_system 
--------------------

As an *alternative* to a ``solvate`` task, you can use a ``make_membrane_system`` task to build a lipid bilayer and *optionally* embed the protein into it.  Pestifer uses ``packmol`` to build a small bilayer "patch" of specified composition, and then replicates it in x and y to generate a sufficiently large bilayer "quilt".

``Pestifer`` comes prepackaged with a set of lipids you can use to build a bilayer.  These lipids are defined in the ``pestifer`` package, and you can see the full list of available lipids by running:

.. code-block:: console

  $ pestifer show-resources --charmmff pdb

``Pestifer`` provides 10 distinct conformers for each lipid molecule which are sampled from short vacuum MD simulations.  These are labelled "0" to "9".  By default, conformer "0" is used for each lipid.

All ``make_membrane_system`` options are documented in the Config Reference pages for :ref:`config_ref tasks make_membrane_system`.

Subtasks
========

There are two main subtask directives in a ``make_membrane_system`` task: ``bilayer`` and ``embed``.  The ``bilayer`` subtask is used to specify the composition of the bilayer, while the (optional) ``embed`` subtask is used to specify how the protein should be embedded in the bilayer.

bilayer
+++++++

This subtask minimally specifies the desired *composition* of the bilayer.  This can be done either using a single ``composition`` directive or by setting packmol-memgen-style string values to the ``lipids`` and ``mole_fractions`` keys.
   
   - **Composition-specification options**

     - ``composition``: Expects a dictionary with two keys: ``lower_leaflet`` and ``upper_leaflet``.  Each of these keys should map to a list of dictionaries, where each dictionary specifies a lipid name and its fraction in the leaflet.  For example, to specify a symmetric bilayer with 50% POPE and 50% POPC in each leaflet, you would use:

       .. code-block:: yaml

         composition:
           lower_leaflet:
             - name: POPE
               frac: 0.5
             - name: POPC
               frac: 0.5
           upper_leaflet:
             - name: POPE
               frac: 0.5
             - name: POPC
               frac: 0.5

       Naturally, you can specify different compositions for the upper and lower leaflets to create an asymmetric bilayer.  The fractions in each leaflet must sum to 1.0.

       You can also specify a variety of conformers for each lipid.  For instance, if you want to use 50% conformer "3" and 50% conformer "5" for POPC, you would specify the composition as follows:

       .. code-block:: yaml

         composition:
           lower_leaflet:
             - name: POPC
               frac: 0.5
               conf: 3
             - name: POPC
               frac: 0.5
               conf: 5
           upper_leaflet:
             - name: POPC
               frac: 0.5
               conf: 3
             - name: POPC
               frac: 0.5
               conf: 5

       Using a variety of conformers is useful for generating a more realistic bilayer structure, as it allows for different lipid conformations to be sampled during the packing process.

      - ``lipids`` and ``mole_fractions``: These strings specify the lipid composition in the bilayer.  Each is a packmol-memgen-style string, where each lipid (or mole fraction) is separated by a comma.  For example, to specify a symmetric bilayer with POPE and POPC, you would use:

        .. code-block:: yaml

          lipids: POPE:POPC
          mole_fractions: 0.5:0.5

        This notation specifies that *both* leaflets have the *same* composition.  For an asymmetric bilayer, the packmole-memgen string format is a bit more complex -- you have to use ``//`` to delineate the two leaflets.  You can specify the upper and lower leaflets separately, like this:

        .. code-block:: yaml

          lipids: POPE:POPC//POPE:POPC
          mole_fractions: 0.5:0.5//0.6:0.4 
        
        You can also use a single packmol-memgen-style string to specify the conformers used for each lipid.  For example, to specify that the bilayer uses conformer "3" and the lower leaflet uses conformer "5" for POPC, you would use:

        .. code-block:: yaml

          lipids: POPC:POPC
          mole_fractions: 0.5:0.5
          conformers: 3:5

   - **Optional subdirectives**: The ``bilayer`` subdirective can also include optional subdirectives to control how the bilayer is built and relaxed:
     
     - ``SAPL`` specifies the target surface area per lipid (SAPL) for the bilayer.  This is a single value that applies to both leaflets.  Pestifer will use this value to determine how many lipids to include in each leaflet of the bilayer patch when using ``packmol``.  Note that this is only an *initial* value; the actual SAPL will be determined during the relaxation protocol specified in the ``relaxation_protocols`` directive.  If you do not specify a ``SAPL``, pestifer will use a default value of 50.0 Angstroms^2 per lipid.
   
     - ``relaxation_protocols`` specifies protocols to relax the bilayer patch and quilt.  Each protocol is a list of tasks, and each task is a dictionary with a single key-value pair.  The key is the name of the task, and the value is a dictionary of options for that task.  The example above shows how to use ``md`` tasks to perform energy minimization, NVT, and NPT equilibration.  The ``patch`` protocol is used to relax the initial bilayer patch, and the ``quilt`` protocol is used to relax the final quilted bilayer system. (``md`` tasks are currently the only supported task type in a ``relaxation_protocols`` directive.  All ``md`` task options are documented in the Config Reference pages for :ref:`config_ref tasks md`.)
     - ``solvents``, ``solvent_mole_fractions``, and  are optional packmol-memgen-style strings that specify the solvent composition in the bilayer patch.  Pestifer uses 100% TIP3P water by default.
     - ``solvent_to_lipid_ratio`` specifies the number of waters to include in the extra-membrane region of the box as proportional to the number of lipids.  If you do not specify this, pestifer will use a default of 32 TIP3P water molecules per lipid in the bilayer patch.
     - ``patch_nlipids`` specifies the number of lipids to include in each bilayer patch leaflet.  It expects a dictionary with keys ``upper`` and ``lower``, each with a value that is an integer.  If you do not specify this, pestifer will use a default of 100 lipids per leaflet in the bilayer patch.  This value is used to determine the size of the bilayer patch when using ``packmol``.
  
``embed``
+++++++++

This subtask requires a set of specifications for embedding the protein in the bilayer.  
   
   - ``xydist`` and ``zdist`` are the lateral and membrane-nomral margins of the simulation box; these values are added to corresponding coordinates extremal atoms to make sure the box is big enough.
   - ``z_head_group`` and ``z_tail_group`` are VMD ``atomselect`` strings that define the center of mass z-coordinate for the membrane-proximal regions of the protein
   - ``z_ref_group`` specifies via VMD atomselection whose center of mass sits at the the middle of the bilayer at a particular ``z_value``.
   - ``no_orient`` is a boolean that specifies whether the protein should be oriented to the membrane normal.  If ``True``, the protein will be oriented so that its z-axis is aligned with the membrane normal.  If ``False``, the protein will not be oriented, and the orientation it has in its own structure file will be preserved.  It is ``False`` by default.  When ``False``, the ``z_head_group`` and ``z_tail_group`` selections are required and used to determine the orientation of the protein in the bilayer.


Example
========

An example ``make_membrane_system`` task is specified below:

.. code-block:: yaml

  - make_membrane_system:
      bilayer:
        SAPL: 50.0
        composition:
          lower_leaflet:
            - name: POPE
              frac: 0.09
            - name: POPE
              frac: 0.09
            - name: SOPS
              frac: 0.18
            - name: SOPE
              frac: 0.30
            - name: CHL1
              frac: 0.43
          upper_leaflet:
            - name: PSM
              frac: 0.36
            - name: POPC
              frac: 0.17
            - name: CHL1
              frac: 0.47
        relaxation_protocols:
          patch:
            md:
              ensemble: minimize
              nsteps: 1000
            md:
              ensemble: NVT
              nsteps: 1000
            md:
              ensemble: NPT
              nsteps: 10000   
          quilt:
            md:
              ensemble: minimize
              nsteps: 1000
            md:
              ensemble: NVT
              nsteps: 1000
            md:
              ensemble: NPT
              nsteps: 1000       
      embed:
        xydist: 30
        zdist: 20
        z_head_group: "protein and resid 667"
        z_tail_group: "protein and resid 710"
        z_ref_group: 
          text: "protein and resid 696"
          z_value: 0.0

If you are including cholesterol or any other sterols in your bilayer, it is recommended that you follow the ``make_membrane_system`` task immediately with an energy minimization and then a ``ring_check`` task.  This will identify and delete any sterols whose rings are pierced by other molecules.  This is illustrated in :ref:`example 14`.

Task Flow
=========

As mentioned above, pestifer first uses packmol to make a minimal patch (of, say, 100 lipids per leaflet) of the desired composition.  If the system is to have a symmetric bilayer (same composition in each leaflet), then this patch is first relaxed and then replicated to form the final quilt, into which the protein is embedded (if there is one).  If the system is to have an asymmetric bilayer (different composition in each leaflet), then pestifer first makes *two symmetric* patches, where in the first the two leaflets have the same composition as the upper leaflet of the final system, and the second has the composition of the lower leaflet.  These two patches are relaxed independently.  Then a hybrid asymmetric patch is constructed by combining the upper leaflet of the first and the lower leaflet of the second.  The lateral box size is set as that of the larger of the two leaflets (laterally).  If there is a difference in lateral area of the two patches, this means the larger one has *excess lipids*.  However, we will not delete excess lipids until the quilt is made.  At this point, this fresh asymmetric patch is replicated to form the quilt.  The number of excess lipids in the larger patch is computed assuming that the equilibrated symmetric patch reports an accurate SAPL for that composition, and that the two leaflets in the quilt must have the *same* area, which is assumed at the outset to reflect a laterally equilibrated *smaller* leaflet.  The determined number of excess lipids is then deleted from the larger leaflet by random selection, and the system is relaxed again.  This is done to ensure the lateral pressure in the two leaflets is the same, minimizing any spontaneous curvature that would cause spurious "ripples" in the fully periodic system.

.. mermaid::
  :caption: Pestifer Make Membrane System Task Bilayer Task Flowchart

  graph TD;
    A{Is bilayer asymmetric?};
    A -- No --> B[Make patch];
    B --> H[Relax];
    H --> I[Replicate patch to quilt];
    A -- Yes --> C[Make two symmetric patches];
    C --> D[Make patch for upper leaflet];
    C --> E[Make patch for lower leaflet];
    D --> F[Relax];
    E --> G[Relax];
    F -- Upper leaflet --> J[Combine leaflets into asymmetric patch];
    G -- Lower leaflet --> J;
    J --> N[Replicate asymmetric patch to quilt];
    N --> K{Any excess lipids?};
    K -- Yes --> L[Delete excess lipids from larger leaflet];
    K -- No --> M[Relax quilted system];
    L --> M;
    I --> M;

Other Notes
===========

Pestifer's ``make_membrane_system`` task is inspired by the `packmol-memgen package <https://ambermd.org/tutorials/advanced/tutorial38/index.php>`_.  For instance, we borrow ``packmol-memgen`` syntax for specifying composition (optionally; the preferred syntax is to use a ``composition`` dictionary in the yaml input).  However, we do not use any precomputed surface-area per lipid.  Instead, we allow the user to specify a single value for SAPL and then a relaxation protocol to achieve a laterally equilibrated bilayer system prior to any embedding.  Packmol-memgen allows for generation of asymmetric bilayers but does not provide a way to guarantee lack of spontaneous curvature that might result, beyond assuming its pre-computed SAPL's are correct.

A moderately sized bilayer system can take *several hours* to pack.  A typical bilayer patch of 100 lipids per leaflet with 32 waters per lipid will take between 15 and 30 minutes for ``packmol`` to pack; the remainder of the time is relaxation.  The relaxation protocols can take anywhere from a few minutes to several hours, depending on the size of the system and the number of steps specified in each protocol.  See :ref:`example 15` for an example of a protein-embedded, heterogeneous asymmetric bilayer system built using the ``make_membrane_system`` task.