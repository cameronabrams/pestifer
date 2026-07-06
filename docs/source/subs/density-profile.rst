.. _subs_density_profile:

density-profile
---------------

The ``density-profile`` subcommand computes and plots species-resolved mass-density profiles — water, total lipid, protein, and ions — along the bilayer normal (:math:`z`) for a membrane system.  It reads a PSF, a single coordinate frame (a PDB or a NAMD binary ``.coor``), and an XSC cell file; the bilayer midplane is centered at :math:`z=0` and bulk solvent is made continuous across the periodic boundary so the water density is not artificially clipped at the box edges.  This is the tool used to produce the density-profile figures for the membrane examples :ref:`16 <example mper-tm symmetric bilayer>` and :ref:`17 <example mper-tm viral bilayer>`.

The simplest invocation names a single ``basename`` and reads ``<basename>.psf``, ``<basename>.coor`` (or ``<basename>.pdb``), and ``<basename>.xsc``:

.. code-block:: bash

   $ pestifer density-profile --basename my_membrane

The files can also be given individually:

.. code-block:: bash

   $ pestifer density-profile --psf sys.psf --coor sys.coor --xsc sys.xsc --out density.png

For a multicomponent bilayer, ``--lipid-components`` decomposes the total lipid density into one curve per lipid species (in addition to the total lipid curve), which exposes leaflet asymmetry — the outer-leaflet species peak on one side, the inner-leaflet species on the other, while cholesterol populates both:

.. code-block:: bash

   $ pestifer density-profile --basename my_membrane --lipid-components

Options
~~~~~~~

- ``--basename`` — basename for ``<basename>.psf``/``.coor`` (or ``.pdb``)/``.xsc`` and, by default, the output image.
- ``--psf`` / ``--coor`` / ``--xsc`` — input files, overriding the basename-derived names.  ``--coor`` accepts either a PDB or a NAMD binary ``.coor``.
- ``--out`` — output PNG (default ``<basename>-density-profile.png``).
- ``--title`` — plot title.
- ``--dz`` — :math:`z`-slab thickness in Å (default ``1.0``).
- ``--lipid-components`` — also plot a profile for each lipid species.
- ``--figsize`` — figure size in inches (default ``6.4 4.4``).
