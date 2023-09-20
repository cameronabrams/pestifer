Welcome to pestifer's documentation!
====================================

**Pestifer** is a Python package for facilitating the use of the `VMD <https://www.ks.uiuc.edu/Research/vmd/>`_ tool `psfgen <https://www.ks.uiuc.edu/Research/vmd/plugins/psfgen/ug.pdf>`_.  Pestifer is intended to make generation of initial coordinates for all-atom MD simulations of large proteins easier.  Yes, CHARMM-GUI is nice, but you might not want to upload your coordinates there, right?

Pestifer development is supported in part by Grants GM100472, AI154071, and AI178833 from the NIH

You provide pestifer a short configuration file and it will retrieve the correct coordinate files from the RCSB and build a solvated, relaxed system.

Pestifer includes the latest CHARMM force field files from the `MacKerell Lab <https://mackerell.umaryland.edu/charmm_ff.shtml>`_ (July 2022).

.. note::

   Pestifer is under active development.

Contents
--------

.. toctree::

   installation
   usage

