.. _config_ref charmmff:

``charmmff``
============

The ``charmff`` directive informs pestifer which CHARMM force-field parameter and topology files are to be used for all psfgen/namd2 executions.  ``standard`` refers to files that are found in the official CHARMMFF release, and are named relative to the topmost ``toppar/`` directory of that release.  ``custom`` refers to files that are not part of the official CHARMMFF release.

Subdirectives:

.. toctree::
   :maxdepth: 1

   charmmff/standard
   charmmff/custom


