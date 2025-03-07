.. _subs_runtasks_bilayer:

bilayer 
-------

As an *alternative* to a ``solvate`` task, you can use a ``bilayer`` task to build a lipid bilayer and embed the protein into it.  Pestifer uses ``packmol`` to build the bilayer and to solvate with TIP3P water.

An example bilayer task is specified below:

.. code-block:: yaml

  - bilayer:
      embed:
        xydist: 30
        zdist: 20
        protein_radius_scaling: 0.9
        z_head_group: "protein and resid 667"
        z_tail_group: "protein and resid 710"
        z_ref_group: 
          text: "protein and resid 696"
          z_value: 0.0
      SAPL: 50.0
      scale_excluded_volume: 0.25
      lipids: POPE:SOPS:SOPE:CHL1//PSM:POPC:CHL1
      mole_fractions: 0.09:0.18:0.30:0.43//0.36:0.17:0.47


There are five main specs in a bilayer task:

1. ``embed`` is a set of specifications for embedding the protein in the bilayer.  
   * ``xydist`` and ``zdist`` are the lateral and membrane-nomral margins of the simulation box; these values are added to corresponding coordinates extremal atoms to make sure the box is big enough.
   * ``z_head_group`` and ``z_tail_group`` are VMD ``atomselect`` strings that define the center of mass z-coordinate for the membrane-proximal regions of the protein, and ``z_ref_group`` specifies via VMD atomselection the middle of the bilayer at a particular z value.
2. ``SAPL`` is the initial surface area per lipid in Å².  This is just used to provide an initial lateral surface area into which the lipid molecules are packed.  
3. ``scale_excluded_volume`` is used to make the lipids "softer" during packing; I haven't really tested whether this is necessary.
4. ``lipids`` allows one to specify the identities of the lipids molecules in each leaflet.  The two leaflets are separated by "//", and each leaflet is represented by a colon-delimited list of lipid molecules, each referenced by their charmff resnames.  We provide a set of PDB files for all lipid molecules defined in the base charmmff topology file ``top_all36_lipid.rtf``.  To see the full list:

.. code-block:: console

   $ pestifer show-resources --charmmff pdb

5. ``mole_fractions`` is where you specify the mole fractions of each molecule in the same order as in the ``lipids`` directive.

A moderately sized bilayer system can take *several hours* to pack.  This is a limitation of using ``packmol``, which currently has no parallel or GPU enabled version.

All ``bilayer`` options are documented in the Config Reference pages for :ref:`config_ref tasks bilayer`.

If you are including cholesterol or any other sterols in your bilayer, it is recommended that you follow the bilayer task immediately with an energy minimization and then a ``ring_check`` task.  This will identify and delete any sterols whose rings are pierced by other molecules.  This is illustrated in :ref:`example 14`.


Pestifer's ``bilayer`` task is superficially inspired by the `packmol-memgen package <https://ambermd.org/tutorials/advanced/tutorial38/index.php>`_. 
