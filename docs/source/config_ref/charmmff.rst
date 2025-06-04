.. _config_ref charmmff:

``charmmff``
============

The ``charmff`` directive informs pestifer which CHARMM force-field parameter and topology files are to be used for all psfgen/namd3 executions.  ``standard`` refers to files that are found in the official CHARMMFF release, and are referred to by their base filenames, regardless of whether they live in the topmost ``toppar`` directory or one of the stream subdirectories.  ``custom`` refers to files that are not part of the official CHARMMFF release.

Subdirectives:

.. toctree::
   :maxdepth: 1

   charmmff/standard
   charmmff/custom
   charmmff/overrides


