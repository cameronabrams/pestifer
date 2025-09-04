.. _introduction:

Introduction
============

Why Pestifer?
-------------

Perhaps the most important phase of large-scale all-atom molecular dynamics (MD) simulations is **system preparation**--starting with an existing atomic structure file and ending with a simulation box and associated topology that can be understood by MD simulation software.  I am a long-time user of `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_ and `NAMD <https://www.ks.uiuc.edu/Research/namd/>`_, and I've invested a lot of time and effort mastering the `psfgen <https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf>`_ VMD plugin to prepare systems that use the `CHARMM36 force field <https://mackerell.umaryland.edu/charmm_ff.shtml>`_.  I have written hundreds of ``psfgen`` scripts over the years, and every new one I wrote made me wish I didn't have to keep reinventing the wheel.  I wanted a way to automate the tedious parts of writing ``psfgen`` scripts, and I wanted a way to ensure that I could reproduce any system I prepared in the future.  So I wrote Pestifer.

Pestifer provides an easy way to generate lots of system replicas using *exactly the same method* each time, and it automates a lot of the tedium encountered in writing scripts for ``psfgen``.  Pestifer requires only a `YAML-format <https://yaml.org/>`_ input configuration file to describe how to prepare a simulation-ready system from an existing structure, in the form of a local structure file (pdb or mmcif formats) or database entry (e.g., RCSB ID or Alphafold ID).  This means that the YAML input together with the identified version of pestifer is a **complete record** of how any system was prepared, and it can be used to reproduce the system at any time in the future.

In addition to enhancing reproducibilty in system preparation, Pestifer also provides for data security.  You never need to upload any data to a web-based application, like `Charmm-GUI <https://charmm-gui.org>`_, to prepare a system.  Pestifer runs entirely on your local machine, so you can prepare simulation systems based on non-publicly-available PDB structures or other data that you do not (yet) want to share with the world.

Why YAML?
---------
`YAML <https://yaml.org/>`_ is a human-readable data format.  Like JSON, it maps neatly to Python containers, but 
unlike JSON, it has support for comments, making it an ideal format for input configuration files for simulations.  It also has a more flexible structure, allowing for complex data representations without the need for extensive punctuation.  I like YAML so much I wrote `Ycleptic <https://ycleptic.readthedocs.io/en/latest/>`_, a Python package that provides a YAML-based configuration system for Python applications.  Pestifer uses Ycleptic to parse its YAML input files.

A partial list of things Pestifer can do
----------------------------------------

1. Generate a complete, self-contained set of files for NAMD simulation of a solvated system using the CHARMM36FF from nearly any PDB input, including glycosylated proteins and nucleic acids
2. During this system-preparation process, introduce modifications, including
   
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

Check out the :ref:`examples` -- some of these capabilities are showcased there.

Who is Pestifer *not* for?
----------------------------

If you are an experienced CHARMM user, then Pestifer is likely not much use to you, since you think in terms of writing a CHARMM input script (or if you hate writing CHARMM scripts, you use ``pyCHARMM``).   If you are unlikely or just unwilling to use NAMD or the CHARMM36 force field, then Pestifer is not for you.  If you are a Charmm-GUI user, then Pestifer may not be for you, since Charmm-GUI has many system-building features that Pestifer does not (yet) have.

Why "Pestifer"?
---------------

The name "Pestifer" is an easy-to-say three-syllable word that includes the letters "P", "S", and "F" exactly once and in that order; since I wanted a reliable way to generate PSF files, I thought that was clever. "Pestifer" comes from the Latin word for "plague-bearing"; you may also adopt the idea that Pestifer is a tool that can help you spread your MD simulations far and wide, like a plague.  Coincidentally, `"Pestifer" <http://pestifer.be/>`_ is also a Belgian technical death-metal band formed in 2004.