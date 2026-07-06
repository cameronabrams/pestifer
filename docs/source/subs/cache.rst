.. _subs_cache:

cache
-----

Pestifer caches expensive-to-parse data so that repeated builds do not have to re-read it each time: the parsed CHARMM force field (topology/parameter data), the built-in PDB repository, and a compact residue-name lookup index.  The caches live in a per-user directory and are keyed to the force-field release; they refresh automatically when the underlying resource files change.  The ``cache`` subcommand lets you inspect and manage them.

status
======

List the caches currently on disk with their sizes and modification times:

.. code-block:: console

    $ pestifer cache status
    pestifer cache directory: /home/you/.cache/pestifer
      charmmffcontent              12.8 MB  2026-07-01 09:13
      charmmffresitopcollection    18.6 MB  2026-07-01 09:51
      pdbrepository                 3.1 MB  2026-07-01 09:50
      resnameindex                 19.2 KB  2026-07-06 11:24
      4 file(s), 34.5 MB total

clear
=====

Delete all cache files.  This is safe -- each cache is rebuilt automatically the next time it is needed (the first build after clearing is slower):

.. code-block:: console

    $ pestifer cache clear

rebuild
=======

Force-rebuild every cache from the current resource files, for each installed CHARMM force-field release.  Use this after changing the packaged force-field files (it is otherwise unnecessary, since the caches refresh on their own when the resources change):

.. code-block:: console

    $ pestifer cache rebuild
