.. _subs_desolvate:

desolvate
---------

`` pestifer desolvate`` generates DCD files that are stripped of solvent (usually water).

``pestifer desolvate`` requires the name of a PSF file and any number of congruent DCD files (listed in chronological order).  It assumes you want to strip out anything that is not protein or lipid, unless you specify a different VMD atomselect string to define the part of ths system you want to keep.  It will generate a new PSF file and a single DCD that is the stripped version of concatenation of the input DCD files.

``pestifer desolvate`` is really just a convenient wrapper for the ``catdcd`` command that comes with VMD.

For help with ``pestifer desolvate``: 

.. code-block:: bash

    $ pestifer desolvate --help

Example
+++++++

Say you have a PSF file ``my_system.psf`` and DCD's generated from sequential NAMD runs called ``run01.dcd``, ``run02.dcd``, etc.  A simple invocation of ``pestifer desolvate`` might be:

.. code-block:: bash

    $ pestifer desolvate --psf my_system.psf --dcd-infiles run??.dcd

This will generate ``dry.psf`` and ``dry.dcd``, which will hold all frames from the input DCD files stripped of everything but protein, lipid and nucleic acids, and any ligands.




