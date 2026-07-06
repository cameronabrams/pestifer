.. _subs_rebuild_charmmff_cache:

rebuild-charmmff-cache
----------------------

Pestifer caches parsed CHARMM force-field topology and parameter data so that repeated builds do not have to re-read and re-parse the full force field each time.  The ``rebuild-charmmff-cache`` subcommand regenerates that cache from the current topology and parameter files.

Run it after adding or updating CHARMM force-field files (for example after installing a new force-field release or adding a custom stream file) so the changes are picked up on the next build:

.. code-block:: bash

   $ pestifer rebuild-charmmff-cache
