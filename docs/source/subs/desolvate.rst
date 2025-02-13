desolvate
---------

Pestifer provides a subcommand called ``desolvate`` to generate DCD files that are stripped of solvent (usually water).

``pestifer desolvate`` requires the name of a PSF file and any number of congruent DCD files (which you will want to be sure are listed in chronological order).  It assumes you want to strip out anything that is not protein or lipid, unless you specify a different VMD atomselect string to define the part of ths system you want to keep.  It will generate a new PSF file and a single DCD that is the stripped version of concatenation of the input DCD files.

``pestifer desolvate`` is really just a convenient wrapper for the ``catdcd`` command.

Example: Say you have a PSF file ``my_system.psf`` and DCD's generated from sequential NAMD runs called ``run01.dcd``, ``run02.dcd``, etc.  A simple invocatino of ``pestifer desolvate`` might be:

.. code-block:: console

    $ pestifer desolvate --psf my_system.psf --dcd-infiles run??.dcd

This will generate ``dry.psf`` and ``dry.dcd``, which will hold all frames from the input DCD files stripped of everything but protein and lipid.

For more help with ``pestifer desolvate``: 

.. code-block:: console

    $ pestifer desolvate --help


