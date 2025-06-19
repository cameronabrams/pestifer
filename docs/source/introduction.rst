.. _introduction:

Introduction
============

Why pestifer?
-------------

``Pestifer`` grew out of my own research in MD simulations of biomolecular systems.  I am a long-time user of `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_ and `NAMD <https://www.ks.uiuc.edu/Research/namd/>`_, and I've invested a lot of time and effort mastering the ``psfgen`` utility to generate NAMD input files for systems that use the CHARMM force field.  I suspect that ``psfgen`` was developed because the CHARMM program was not freely available but the CHARMM force field was; this is likely why there are so many of us who think ``psfgen/VMD`` when building systems that use the CHARMM force-field and are run using NAMD.  

If you are an experienced CHARMM user, then ``pestifer`` is likely not much use to you, since you think in terms of writing a CHARMM input script (or if you hate writing CHARMM scripts, you use ``pyCHARMM``). If you are a regular user of `Charmm-GUI <https://charmm-gui.org>`_, a server-based interface to the CHARMM PSF-generation routines that can generate production-ready NAMD topology and coordinate inputs, ``pestifer`` offers some advantages.  Unless you are exclusively using coordinates deposited at the `RCSB Protein Data Bank <https://rcsb.org>`_, Charmm-GUI requires you to upload a PDB file to begin the process, and maybe you aren't comfortable with that.  ``pestifer`` does everything on your local machine, so you can build simulation systems based on non-publicly-available PDB structures.  ``pestifer`` also works well in a SLURM-based HPC environment.

I wrote pestifer to provide an easy way to generate lots of system replicas using *exactly the same method* each time, and to automate a lot of the tedium encountered in writing scripts for ``psfgen``.  Pestifer requires only a YAML-format input configuration file to describe how to build a simulation-ready system from an existing structure, in the form of a local structure file (pdb or mmcif) or database entry (e.g., RCSB PDB ID or Alphafold ID).  This means that the YAML input together with the identified version of pestifer is a complete record of how any system was built, and it can be used to reproduce the system at any time in the future.  This is the key reason for pestifer.

As far as automating the tedium of psfgen scripting, pestifer can do a lot of things:

A partial list of things pestifer can do 
----------------------------------------

1. Generate a complete, self-contained set of files for NAMD simulation of a solvated system using the CHARMM36FF from nearly any PDB input
2. During this generation process, introduce modifications, including
   
   * residue patches
   * residue mutations
   * residue insertions
   * residue deletions
   * residue substitutions
   * fusions
   * automatic construction of biomolecular assemblies
   * automatic loop-building to model unresolved residues
   * arbitrary backbone and side-chain rotations
   * grafting of pendant groups from one PDB onto another (for glycans, mostly)
   * chain cleavage
   * domain-swapping
   * chain relabeling
3. Generate membrane-embedded proteins using any lipids for which a PDB is available (using packmol), including automatic checking for pierced rings
4. Generate restart files (NAMD configs and SLURM scripts) based on current runs 
5. Generate solvent-stripped PSF and DCD files

Check out the :ref:`examples` -- some of these are showcased there.
