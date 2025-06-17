.. _subs_runtasks_psfgen:

psfgen 
------

A ``psfgen`` task is usually the first task in a ``tasks`` directive.  Its basic functionality is to set up and conduct the first ``psfgen`` run to generate psf and pdb files from a structure file.  Because pestifer imports ``pidibble``, it can also fetch structure files from the RSCB or AlphaFold.

An example ``psfgen`` task that fetches and processes the 6PTI structure might look like this:

.. code-block:: yaml

   psfgen:
    source:
      id: 6pti


With this input, pestifer will fetch ``6pti.pdb`` from the RCSB, generate the required ``psfgen`` script, and execute ``psfgen`` from within a text-only VMD session.

The only two subdirectives under ``psfgen`` are :ref:`config_ref tasks psfgen source` and :ref:`config_ref tasks psfgen mods`.

Under ``source``, there are three mutually exclusive directives used to specify the molecular input:

1. ``id``: This is a 4-character PDB ID.
2. ``alphafold``: This is an AlphaFold database ID.
3. :ref:`config_ref tasks psfgen source prebuilt`: This allows you to specify an existing set containing a psf, pdb, and xsc file for a system that was already built.

Secondary directives under ``source`` are described in the Config File reference pages at :ref:`config_ref tasks psfgen source`. 

Under ``mods``, you can specify any number of modifications that can be handled in a single invocation of ``psfgen``.  These include point mutations, deletions, substitutions, additional disulfide bonds or other covalent bonds, among others.  

For example, a ``psfgen`` task that builds a 6PTI system with two specific point mutations might look like this:

.. code-block:: yaml

   psfgen:
     source:
       id: 6pti
     mods:
       mutations:
         - A:T11A # threonine to alanine at position 11, 
         - A:PRO,13,ALA # alternate syntax; proline to alanine at position 13

Below are all the ``mods`` that can be invoked in a ``mods`` directive.  Each of these is described in detail in the Config File reference pages at :ref:`config_ref tasks psfgen mods`.

.. toctree::
   :maxdepth: 1

   mods/patches
   mods/mutations
   mods/deletions
   mods/insertions
   mods/crotations
   mods/transrot
   mods/ssbonds
   mods/links
   mods/substitutions
   mods/grafts
