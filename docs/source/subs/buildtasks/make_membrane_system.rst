.. _subs_buildtasks_make_membrane_system:

make_membrane_system 
--------------------

As an *alternative* to a ``solvate`` task, you can use a ``make_membrane_system`` task to build a lipid bilayer and *optionally* embed the protein into it.

Pestifer generates the initial lipid and solvent coordinates with a grid packer: it places lipids on a per-leaflet 2-D lattice (oriented, head groups out) and solvent on a 3-D lattice **directly**, and relies on the downstream relaxation MD to resolve the initial overlaps.  It builds the full-size membrane in one shot -- sized to the embedded protein's xy footprint plus the ``embed.xydist`` margin -- so there is no separate patch-packing or replication step.  The `Task Flow`_ section below describes the process.

``Pestifer`` comes prepackaged with a set of lipid PDB files you can use to build a bilayer.  You can see the full list of available lipids witht this command:

.. code-block:: console

  $ pestifer show-resources pdb-repo

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
     
     - ``SAPL`` specifies the target surface area per lipid (SAPL) for the bilayer, applied to both leaflets.  It sets the initial lattice spacing.  This is only an *initial* value; the equilibrium area is found during the ``relaxation_protocols`` relaxation.  If you do not specify a ``SAPL``, pestifer will use a default value of 75.0 :math:`Å^2` per lipid.
   
     - ``relaxation_protocols`` specifies protocols to relax the bilayer before embedding.  Each protocol is a list of ``md`` tasks (the only supported task type), each a dictionary with a single ``md`` key whose value is a dictionary of options.  The ``patch`` protocol relaxes the two per-leaflet *calibration* patches (asymmetric builds only), and the ``quilt`` protocol relaxes the full-size, directly-gridded membrane.  Long ``NPT``/``NPgT`` stages are automatically split into a ramp of restarts so the condensing cell never outgrows NAMD's startup patch grid.  All ``md`` task options are documented in the Config Reference pages for :ref:`config_ref tasks md`.
     - ``solvents`` and ``solvent_mole_fractions`` are optional packmol-memgen-style strings that specify the solvent composition.  Pestifer uses 100% TIP3P water by default.
     - ``solvent_to_lipid_ratio`` specifies the number of waters to include in the extra-membrane region of the box as proportional to the number of lipids.  If you do not specify this, pestifer will use a default of 32 TIP3P water molecules per lipid.
     - ``patch_nlipids`` specifies the number of lipids in each calibration-patch leaflet (asymmetric builds only).  It expects a dictionary with keys ``upper`` and ``lower``, each an integer.  If you do not specify this, pestifer will use a default of 100 lipids per leaflet.
  
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

If you are including cholesterol or any other sterols in your bilayer, a ``make_membrane_system`` build now guards against pierced sterol rings automatically: during the membrane's relaxation, just before the first dynamics stage, it inserts a lipid ``ring_check`` (which deletes the offending lipid) followed by a short minimize, so a lipid tail threaded through a sterol ring is removed before it can destabilize the run.  It is still recommended to follow the ``make_membrane_system`` task with an energy minimization and an explicit ``ring_check`` to catch glycan or protein rings (which the automatic guard, checking only lipids, does not) when the embedded protein is glycosylated.  This is illustrated in :ref:`example mper-tm viral bilayer`.

Task Flow
=========

.. note::

   "Quilt" is a legacy label for the *full-size* membrane and its relaxation -- it
   dates from an earlier packer that tiled small patches into a quilt.  The grid
   packer does no tiling, so a build's ``*-quilt`` artifacts (``grid-quilt``,
   ``psfgen-quilt``, ``equilibration-quilt``) are the whole membrane, **not**
   patches.  Patches appear only in *asymmetric* builds, as the per-leaflet
   calibration patches described below.

For a **symmetric** bilayer, pestifer builds the full-size membrane directly: it grids both leaflets at the requested composition and initial ``SAPL``, sizing the box to the embedded protein's xy footprint plus the ``embed.xydist`` margin (or to ``patch_nlipids`` otherwise).  The gridded membrane is relaxed with the ``quilt`` protocol, and the protein is then embedded.  There is no separate patch or replication step.

For an **asymmetric** bilayer (different leaflet compositions), pestifer first relaxes two symmetric *calibration* patches -- one per leaflet composition -- with the ``patch`` protocol, to measure each leaflet's preferred area per lipid (APL).  It then grids the full membrane (sized to the protein footprint when embedding) at per-leaflet counts in the stress-free ratio :math:`n_\text{upper}/n_\text{lower} = \text{APL}_\text{lower}/\text{APL}_\text{upper}`, so the two leaflets carry equal area at zero differential stress *by construction* -- no leaflet extraction or excess-lipid deletion is required.  The full membrane is then relaxed with the ``quilt`` protocol before embedding.

.. mermaid::
  :caption: Membrane build flow.

  graph TD;
    A{Is bilayer asymmetric?};
    A -- No --> B[Grid full membrane at initial SAPL];
    A -- Yes --> C[Relax two per-leaflet calibration patches];
    C --> D[Measure each leaflet's preferred APL];
    D --> E[Grid full membrane at stress-free per-leaflet counts];
    B --> F[Relax quilt];
    E --> F;
    F --> G[Embed protein];

Other Notes
===========

Pestifer's ``make_membrane_system`` task is inspired by the `packmol-memgen package <https://ambermd.org/tutorials/advanced/tutorial38/index.php>`_.  For instance, we borrow ``packmol-memgen`` syntax for specifying composition (optionally; the preferred syntax is to use a ``composition`` dictionary in the yaml input).  However, we do not use any precomputed surface-areas per lipid.  Instead, we allow the user to specify a single value for SAPL and then use an MD-based relaxation protocol to achieve a laterally equilibrated bilayer system prior to any embedding.  Packmol-memgen allows for generation of asymmetric bilayers but does not provide a way to guarantee lack of spontaneous curvature that might result, beyond assuming its pre-computed SAPL's are correct.

Most of the wall-clock time for a membrane build is relaxation MD, which can run from a few minutes to several hours depending on the system size and the number of steps in each protocol.  The *coordinate generation* itself is cheap: the grid packer places lipids deterministically in seconds.  See :ref:`example mper-tm viral bilayer` for an example of a protein-embedded, heterogeneous asymmetric bilayer system.