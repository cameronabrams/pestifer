.. _subs_runtasks_make_membrane_system:

make_membrane_system 
--------------------

As an *alternative* to a ``solvate`` task, you can use a ``make_membrane_system`` task to build a lipid bilayer and optionally embed the protein into it.  Pestifer uses ``packmol`` to build a small bilayer "patch" of specified composition and symmetry, and then replicates it in x and y to generate a sufficiently large bilayer "quilt".

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


There are two main directives in a ``make_membrane_system`` task:

1. ``bilayer`` is a set of specifications for the composition of the bilayer that allows packmol to build it
   
   - ``SAPL`` is the initial surface area per lipid in Å².  This is just used to provide an initial lateral surface area into which the lipid molecules are packed.
   - ``composition`` allows one to specify the identities and mole fractions of the lipids molecules in each leaflet.  We provide a set of PDB files for all lipid molecules defined in the base charmmff topology file ``top_all36_lipid.rtf`` along with a few select lipids defined in CHARMM lipid stream files.  To see the full list:

   .. code-block:: console

      $ pestifer show-resources --charmmff pdb

   - ``relaxation_protocols`` is a set of protocols that are used to relax the bilayer patch and quilt.  Each protocol is a list of tasks, and each task is a dictionary with a single key-value pair.  The key is the name of the task, and the value is a dictionary of options for that task.  The example above shows how to use ``md`` tasks to perform energy minimization, NVT, and NPT equilibration.  The ``patch`` protocol is used to relax the initial bilayer patch, and the ``quilt`` protocol is used to relax the final quilted bilayer system. (``md`` tasks are currently the only supported task type in a ``relaxation_protocols`` directive.  All ``md`` task options are document in the Config Reference pages for :ref:`config_ref tasks md`.)

2. ``embed`` is a set of specifications for embedding the protein in the bilayer.  
   
   - ``xydist`` and ``zdist`` are the lateral and membrane-nomral margins of the simulation box; these values are added to corresponding coordinates extremal atoms to make sure the box is big enough.
   - ``z_head_group`` and ``z_tail_group`` are VMD ``atomselect`` strings that define the center of mass z-coordinate for the membrane-proximal regions of the protein, and ``z_ref_group`` specifies via VMD atomselection the middle of the bilayer at a particular z value.

A moderately sized bilayer system can take *several hours* to pack.  This is a limitation of using ``packmol``, which currently has no parallel or GPU enabled version.

All ``make_membrane_system`` options are documented in the Config Reference pages for :ref:`config_ref tasks make_membrane_system`.

If you are including cholesterol or any other sterols in your bilayer, it is recommended that you follow the ``make_membrane_system`` task immediately with an energy minimization and then a ``ring_check`` task.  This will identify and delete any sterols whose rings are pierced by other molecules.  This is illustrated in :ref:`example 14`.

As mentioned above, pestifer first uses packmol to make a minimal patch (of, say, 100 lipids per leaflet) of the desired composition.  If the system is to have a symmetric bilayer (same composition in each leaflet), then this patch is first relaxed and then replicated to form the final quilt, into which the protein is embedded (if there is one).  If the system is to have an asymmetric bilayer (different composition in each leaflet), then pestifer first makes *two symmetric* patches, where in the first the two leaflets have the same composition as the upper leaflet of the final system, and the second has the composition of the lower leaflet.  These two patches are relaxed.  Then a hybrid asymmetric patch is constructed by combining the upper leaflet of the first and the lower leaflet of the second.  The lateral box size is set as that of the larger of the two leaflets (laterally).  If there is a difference in lateral area of the two patches, this means the larger one has *excess lipids*.  However, we will not delete excess lipids until the quilt is made.  At this point, this fresh asymmetric patch is replicated to form the quilt.  The number of excess lipids in the larger patch is computed assuming that the equilibrated symmetric patch reports an accurate SAPL for that composition, and that the two leaflets in the quilt must have the *same* area, which is assumed at the outset to reflect a laterally equilibrated *smaller* leaflet.  The excess lipids are then deleted from the larger leaflet, and the system is relaxed again.  This is done to ensure the lateral pressure in the two leaflets is the same, minimizing any spontaneous curvature that would cause spurious "ripples" in the fully periodic system.

Pestifer's ``make_membrane_system`` task is inspired by the `packmol-memgen package <https://ambermd.org/tutorials/advanced/tutorial38/index.php>`_.  For instance, we borrow ``packmol-memgen`` syntax for specifying composition.  However, we do not use any precomputed surface-area per lipid.  Instead, we allow the user to specify a single value for SAPL and then a relaxation protocol to achieve a laterally equilibrated bilayer system prior to any embedding.  Packmol-memgen allows for generation of asymmetric bilayers but does not provide a way to guarantee lack of spontaneous curvature that might result, beyond assuming its pre-computed SAPL's are correct.
