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

prebuild
========

Generate one lipid conformer set into the cache **before** a build asks for it.

A build that needs a ``(lipid, phase)`` pair with no entry generates that entry itself and caches it under ``~/.pestifer/pdbrepository/<release>/lipid/``, so this is never required.  It is worth doing when generating inside a job is awkward: the generation is a single-molecule vacuum run, so on a cluster it spends allocation time on work that a login node can do, and it happens while the rest of the build waits.

.. code-block:: console

    $ pestifer cache prebuild --resname PSM --phase Lo
    INFO> prebuilding conformer set PSM__Lo (mc sampler)
    ...
    /home/you/.pestifer/pdbrepository/feb26/lipid/PSM__Lo

``--phase`` is ``Ld`` (the fluid ensemble, cached under the bare ``<RESI>`` name) or ``Lo`` (the ordered one, cached as ``<RESI>__Lo``); ``--charmmff-release`` selects the release to build against, defaulting to the newest installed, which is what a build uses.

The sampler is not a choice here.  It is the one a build would pick for that phase, so the cached entry is the entry the build would have made -- a set generated another way would be a cache hit that silently changes what gets packed.

Run it once per lipid and phase your compositions name.  A leaflet like

.. code-block:: yaml

    composition:
      lower_leaflet_phase: Lo
      lower_leaflet:
        - {name: PSM, frac: 0.36}
        - {name: POPC, frac: 0.17}
        - {name: CHL1, frac: 0.47}

resolves to the three entries ``PSM__Lo``, ``POPC__Lo`` and ``CHL1__Lo``.
