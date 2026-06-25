.. _subs_buildtasks_make_membrane_system:

make_membrane_system 
--------------------

As an *alternative* to a ``solvate`` task, you can use a ``make_membrane_system`` task to build a lipid bilayer and *optionally* embed the protein into it.

Pestifer offers two ways to generate the initial lipid and solvent coordinates, selected with the ``bilayer.packer`` option:

- ``grid`` places lipids on a per-leaflet 2-D lattice (oriented, head groups out) and solvent on a 3-D lattice **directly**, and relies on the downstream relaxation MD to resolve the initial overlaps.  It builds the full-size membrane in one shot -- sized to the embedded protein's xy footprint plus the ``embed.xydist`` margin -- so there is no separate patch-packing or replication step, it needs **no** ``packmol`` installation (``paths.packmol`` is not consulted), and it is orders of magnitude faster than packmol.  Examples 16 and 17 use the grid packer.
- ``packmol`` (the default) uses `packmol <https://m3g.github.io/packmol/>`_ to pack a small bilayer "patch" of the specified composition and then replicates it in x and y to form a sufficiently large "quilt".  It requires a ``packmol`` installation (see :ref:`installation <installation>`).

Both packers feed the same relaxation and embedding steps; they differ only in how the starting coordinates are produced.  The `Task Flow`_ section below describes each.

``Pestifer`` comes prepackaged with a set of lipid PDB files you can use to build a bilayer.  You can see the full list of available lipids witht this command:

.. code-block:: console

  $ pestifer show-resources --charmmff pdb

``Pestifer`` provides 10 distinct conformers for each lipid molecule which are sampled from short vacuum MD simulations.  These are labelled "0" to "9".  By default, conformer "0" is used for each lipid.

All ``make_membrane_system`` options are documented in the Config Reference pages for :ref:`config_ref tasks make_membrane_system`.

Subtasks
========

There are two main subtask directives in a ``make_membrane_system`` task: ``bilayer`` and ``embed``.  The ``bilayer`` subtask requires specification of the composition (and optionally the size) of the bilayer, while the (optional) ``embed`` subtask requires specification of how the protein should be embedded in the bilayer.

bilayer
+++++++

Bilayer composition is specified either using a single ``composition`` directive or by setting ``packmol-memgen``-style string values to the ``lipids`` and ``mole_fractions`` keys.
   
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

      - ``lipids`` and ``mole_fractions``: These strings specify the lipid composition in the bilayer.  Each is a ``packmol-memgen``-style string, where each lipid (or mole fraction) is separated by a comma.  For example, to specify a symmetric bilayer with POPE and POPC, you would use:

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
     
     - ``packer`` selects the coordinate generator: ``packmol`` (the default) or ``grid`` (see the introduction above and `Task Flow`_).  The grid packer needs no ``packmol`` installation.
     - ``SAPL`` specifies the target surface area per lipid (SAPL) for the bilayer, applied to both leaflets.  It sets how many lipids fill each leaflet (packmol patch) or the initial lattice spacing (grid).  This is only an *initial* value; the equilibrium area is found during the ``relaxation_protocols`` relaxation.  If you do not specify a ``SAPL``, pestifer will use a default value of 75.0 :math:`Å^2` per lipid.
   
     - ``relaxation_protocols`` specifies protocols to relax the bilayer before embedding.  Each protocol is a list of ``md`` tasks (the only supported task type), each a dictionary with a single ``md`` key whose value is a dictionary of options.  The ``patch`` protocol relaxes the small patches -- the packmol patch, or the grid packer's two per-leaflet *calibration* patches (asymmetric builds only) -- and the ``quilt`` protocol relaxes the full-size membrane (the replicated quilt for packmol, or the directly-gridded membrane for grid).  Long ``NPT``/``NPgT`` stages are automatically split into a ramp of restarts so the condensing cell never outgrows NAMD's startup patch grid.  All ``md`` task options are documented in the Config Reference pages for :ref:`config_ref tasks md`.
     - ``solvents`` and ``solvent_mole_fractions`` are optional packmol-memgen-style strings that specify the solvent composition.  Pestifer uses 100% TIP3P water by default.
     - ``solvent_to_lipid_ratio`` specifies the number of waters to include in the extra-membrane region of the box as proportional to the number of lipids.  If you do not specify this, pestifer will use a default of 32 TIP3P water molecules per lipid.
     - ``patch_nlipids`` specifies the number of lipids in each patch leaflet (the packmol patch, or the grid calibration patches).  It expects a dictionary with keys ``upper`` and ``lower``, each an integer.  If you do not specify this, pestifer will use a default of 100 lipids per leaflet.
  
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

If you are including cholesterol or any other sterols in your bilayer, it is recommended that you follow the ``make_membrane_system`` task immediately with an energy minimization and then a ``ring_check`` task.  This will identify and delete any sterols whose rings are pierced by other molecules.  This is illustrated in :ref:`example mper-tm viral bilayer`.

Task Flow
=========

The detailed flow depends on the ``packer``.

Grid packer
+++++++++++

For a **symmetric** bilayer, the grid packer builds the full-size membrane directly: it grids both leaflets at the requested composition and initial ``SAPL``, sizing the box to the embedded protein's xy footprint plus the ``embed.xydist`` margin (or to ``patch_nlipids`` otherwise).  The gridded membrane is relaxed with the ``quilt`` protocol, and the protein is then embedded.  There is no separate patch or replication step.

For an **asymmetric** bilayer (different leaflet compositions), the grid packer first relaxes two symmetric *calibration* patches -- one per leaflet composition -- with the ``patch`` protocol, to measure each leaflet's preferred area per lipid (APL).  It then grids the full membrane (sized to the protein footprint when embedding) at per-leaflet counts in the stress-free ratio :math:`n_\text{upper}/n_\text{lower} = \text{APL}_\text{lower}/\text{APL}_\text{upper}`, so the two leaflets carry equal area at zero differential stress *by construction* -- no leaflet extraction or excess-lipid deletion is required.  The full membrane is then relaxed with the ``quilt`` protocol before embedding.

.. mermaid::
  :caption: Grid packer flow.

  graph TD;
    A{Is bilayer asymmetric?};
    A -- No --> B[Grid full membrane at initial SAPL];
    A -- Yes --> C[Relax two per-leaflet calibration patches];
    C --> D[Measure each leaflet's preferred APL];
    D --> E[Grid full membrane at stress-free per-leaflet counts];
    B --> F[Relax quilt];
    E --> F;
    F --> G[Embed protein];

Packmol packer
++++++++++++++

The packmol packer first uses packmol to make a minimal patch (of, say, 100 lipids per leaflet) of the desired composition.  If the system is to have a symmetric bilayer (same composition in each leaflet), then this patch is first relaxed and then replicated to form the final quilt, into which the protein is embedded (if there is one).  If the system is to have an asymmetric bilayer (different composition in each leaflet), then pestifer first makes *two symmetric* patches, where in the first the two leaflets have the same composition as the upper leaflet of the final system, and the second has the composition of the lower leaflet.  These two patches are relaxed independently.  Then a hybrid asymmetric patch is constructed by combining the upper leaflet of the first and the lower leaflet of the second.  The lateral box size is set as that of the larger of the two leaflets (laterally).  If there is a difference in lateral area of the two patches, this means the larger one has *excess lipids*.  However, we will not delete excess lipids until the quilt is made.  At this point, this fresh asymmetric patch is replicated to form the quilt.  The number of excess lipids in the larger patch is computed assuming that the equilibrated symmetric patch reports an accurate SAPL for that composition, and that the two leaflets in the quilt must have the *same* area, which is assumed at the outset to reflect a laterally equilibrated *smaller* leaflet.  The determined number of excess lipids is then deleted from the larger leaflet by random selection, and the system is relaxed again.  This is done to ensure the lateral pressure in the two leaflets is the same, minimizing any spontaneous curvature that would cause spurious "ripples" in the fully periodic system.

.. mermaid::
  :caption: Packmol packer flow.

  graph TD;
    A{Is bilayer asymmetric?};
    A -- No --> B[Pack patch];
    B --> H[Relax];
    H --> I[Replicate patch to quilt];
    A -- Yes --> C[Pack two symmetric patches];
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

Pestifer's ``make_membrane_system`` task is inspired by the `packmol-memgen package <https://ambermd.org/tutorials/advanced/tutorial38/index.php>`_.  For instance, we borrow ``packmol-memgen`` syntax for specifying composition (optionally; the preferred syntax is to use a ``composition`` dictionary in the yaml input).  However, we do not use any precomputed surface-areas per lipid.  Instead, we allow the user to specify a single value for SAPL and then use an MD-based relaxation protocol to achieve a laterally equilibrated bilayer system prior to any embedding.  Packmol-memgen allows for generation of asymmetric bilayers but does not provide a way to guarantee lack of spontaneous curvature that might result, beyond assuming its pre-computed SAPL's are correct.

Most of the wall-clock time for a membrane build is relaxation MD, which can run from a few minutes to several hours depending on the system size and the number of steps in each protocol.  The two packers differ markedly in how long the *coordinate generation* itself takes: the ``grid`` packer places lipids deterministically in seconds, whereas a packmol patch of 100 lipids per leaflet with 32 waters per lipid takes ``packmol`` roughly 15-30 minutes to pack.  See :ref:`example mper-tm viral bilayer` for an example of a protein-embedded, heterogeneous asymmetric bilayer system built using the grid packer.