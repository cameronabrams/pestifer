Welcome to pestifer's documentation!
====================================

**Pestifer** is a Python package that facilitates the use of the `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_ tool `psfgen <https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf>`_ to generate CHARMM-force-field compatible PSF and PDB files for use in the MD simulation package NAMD.  Pestifer automates and extends the ``psfgen/VMD`` system-building process using a simple YAML interface.

Pestifer grew out of my own research in MD simulations of biomolecular systems as a user of VMD and `NAMD <https://www.ks.uiuc.edu/Research/namd/>`_, and the need to carefully document how systems are generated.  I suspect that ``psfgen`` was developed because the CHARMM program was not freely available but the CHARMM force-field was; this is likely why there are so many of us who think ``psfgen/VMD`` when building systems that use the CHARMM force-field and are run using NAMD.  If you are an experienced CHARMM user, then pestifer is likely not much use to you since you think in terms of writing a CHARMM input script (or if you hate writing CHARMM scripts, you use ``pyCHARMM``). If you are a regular user of `Charmm-GUI <https://charmm-gui.org>`_, a server-based interface to the CHARMM PSF-generation routines that can generate production-ready NAMD topology and coordinate inputs, pestifer offers some advantages.  Unless you are exclusively using coordinates deposited at the `RCSB Protein Data Bank <https://rcsb.org>`_, Charmm-GUI requires you to upload a PDB file to begin the process, and maybe you aren't comfortable with that.  ``pestifer`` does everything on your local machine, so you can build simulation systems based on non-publicly-available PDB structures.

Pestifer includes the CHARMM force field files from the `MacKerell Lab <https://mackerell.umaryland.edu/charmm_ff.shtml>`_ (July 2022).

.. note::

   Pestifer is under active development.

.. note::

   Pestifer development is supported in part by Grants GM100472, AI154071, and AI178833 from the NIH.

Contents
--------

.. toctree::

   installation
   usage
   examples

