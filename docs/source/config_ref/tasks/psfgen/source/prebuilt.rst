prebuilt 
--------

This directive instructs pestifer to load an existing system which already has minimally a PSF and PDB file, and optionally an XSC file.

* ``psf``: name of PSF file
* ``pdb``: name of PDB file
* ``xsc``: (optional) name of XSC file

Example:

.. code-block:: yaml

    tasks:
      psfgen:
        prebuilt:
          pdb: my_build.pdb
          psf: my_build.psf
          xsc: solvated.xsc
    