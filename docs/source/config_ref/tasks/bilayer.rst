.. _config_ref tasks bilayer:

``bilayer``
===========

Parameters controlling packmol to generate a bilayer system

Single-valued parameters:

  * ``lipids``: declaration of lipids in each leaflet in the format LIP1[:LIP2][//LIP3] (default: POPC)

  * ``lipid_conformers``: declaration of which conformer index to use for each lipid in the format C1:C2//C3 (defaults to 0) (default: 0)

  * ``mole_fractions``: declaration of mole fractions of each lipid in each leaflet in the format R1:R2//R3 (default: 1.0)

  * ``solvents``: string in form of SOL1:SOL2:... (default: TIP3)

  * ``solvent_mole_fractions``: string in form of x_1:x_2:... (default: 1.0)

  * ``SAPL``: estimate of surface area per lipid (A^2) (default: 60.0)

  * ``leaflet_thickness``: initial value of leaflet thickness (A) (default: 20.0)

  * ``scale_excluded_volume``: scaling factor controlling contribution of non-solvent atoms in calculating available chamber volumes (default: 1.0)

  * ``fuzz_factor``: when tilted a little from the global-z axis, a long molecule longest z-dimension is lower than its longest overall dimension; this difference is called "fuzz".  When packmol places a long molecule with some reference atoms above a z-plane and other reference atoms below that z-plane, this factor is used to assign the "fuzz" to either side of the planes (this is a terrible explanation; see docs) (default: 0.5)

  * ``dims``: box dimensions (A); must be specified for a bilayer-only system

  * ``length_pad``: additional length added to longest molecular length when orienting molecules that are parallel to z axis

  * ``solution_gcc``: solution density in g/cc (default: 1.0)

  * ``cation``: name of salt cation (default: POT)

  * ``anion``: name of salt anion (default: CLA)

  * ``salt_con``: salt concentration in M

  * ``nloop``: number of packmol GENCAN loops for every component

  * ``nloop_all``: number of packmol GENCAN loops for altogether packing

  * ``tolerance``: clash detection tolerance for packmol

  * ``seed``: RNG seed (optional)



Subdirective:

.. toctree::
   :maxdepth: 1

   bilayer/embed


