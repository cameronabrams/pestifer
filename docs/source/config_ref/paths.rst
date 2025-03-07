.. _config_ref paths:

``paths``
=========

Allows user to specify paths for certain resources, including external applications

Single-valued parameters:

  * ``namd2``: Path for namd2 executable (default: namd2)

  * ``namd3``: Path for namd3 executable, CPU version (default: namd3)

  * ``namd3gpu``: Path for namd3 executable, GPU version; this must be different from the path for the CPU version if both are to be used in a single Pestifer build (default: namd3gpu)

  * ``charmrun``: Path for charmrun executable (default: charmrun)

  * ``vmd``: Path for vmd executable (default: vmd)

  * ``charmmff``: Path for standard charmm force-field distribution; note that pestifer includes the July 2024 distribution by default

  * ``packmol``: Path to packmol executable (default: packmol)

  * ``catdcd``: Path to catdcd executable (default: catdcd)

  * ``pdb_depot``: Path to a user's collection of PDB files for custom lipids or ligands



