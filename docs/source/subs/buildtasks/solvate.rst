.. _subs_buildtasks_solvate:

solvate 
-------

This task simply performs a standard ``psfgen`` run that invokes the ``solvate`` and ``autoionize`` VMD plugins to make a system solvated by TIP3P water.  It is typically specified using its name only.  A prototypical usage example is below:

.. code-block:: yaml

   tasks:
     - fetch:
         sourceID: 6pti
     - psfgen:
     - md:
         ensemble: minimize
     - md:
         ensemble: NVT
     - solvate:
         salt_con: 0.15
     - md:
          ensemble: minimize
     - md:
         ensemble: NVT
    
After minimizing and thermalizing (300 K) a newly built system, pestifer will then solvate it and autoionize to achieve a salt concentration of 150 mM, and then it performs a second minimization and thermalization.

Solvating with a non-water solvent
++++++++++++++++++++++++++++++++++

By default the ``solvate`` task fills with TIP3P water (VMD's built-in pre-equilibrated box).  To solvate with a different solvent, set the ``solvent`` key to the solvent's RESI name:

.. code-block:: yaml

   - solvate:
       solvent: MEOH
       salt_con: 0.15

Any ``solvent`` other than ``TIP3``/``water`` is looked up as a ``kind: box`` entry in the ``solvent`` PDB collection: pestifer checks out the entry's pre-equilibrated box psf/pdb and hands them to VMD's ``solvate`` plugin as ``-spsf``/``-spdb`` together with the box's exact equilibrated edge (``-ws``) and key atom (``-ks``), so the box tiles the simulation cell instead of the built-in water box.  ``autoionize`` still runs afterward, so ``salt_con``/``cation``/``anion`` behave as usual.

Pestifer ships pre-equilibrated boxes for **methanol (**\ ``MEOH``\ **), ethanol (**\ ``ETOH``\ **), and DMSO (**\ ``DMSO``\ **)**, so those work with no extra setup.

For **any other** solvent that is defined in the CHARMM force field, pestifer **builds a box for it on demand** and caches the result per-user, so ``solvent: <RESI>`` just works:

.. code-block:: yaml

   solvate:
     solvent: BENZ        # not shipped -> pestifer builds + caches a benzene box, then proceeds

The first build is a one-time cost (pack + minimize + NPT equilibration, loudly logged) and the box is cached under ``~/.pestifer/pdbrepository/<release>/solvent/<RESI>/``; every subsequent build reuses it instantly.  The cache is keyed by CHARMMFF release (parameters are release-specific) and auto-registered as a user PDB collection.  This *complements* the contribute flow below — the cache is a per-user convenience; ``make-pdb-collection`` + ``modify-package`` is still how you get a curated box into the *shipped* package.

To disable on-demand generation (for reproducibility-strict or offline runs), set ``charmmff.generate_missing_coordinates: false``; a missing box then hard-errors, and you must **build and install** it yourself (see :ref:`sub_make_pdb_collection`):

.. code-block:: bash

   $ pestifer make-pdb-collection solvent --resname MEOH --nmol 216 --density 0.792
   $ pestifer modify-package pdb-repo add-entry solvent/MEOH --collection solvent

.. note::

   Building a solvent box needs a 3-D copy of the solvent molecule.  Water and ions have PDB-repository entries and most residues carry an internal-coordinate table; for small-molecule CGenFF solvents that have neither (e.g. ``MEOH``, ``ETOH``, ``DMSO``), pestifer falls back to **Open Babel** (``obabel --gen3d``) to generate coordinates from the RESI's atom+bond graph.  This requires ``obabel`` on the ``PATH``.  See :ref:`sub_make_pdb_collection`.

Ions and neutralization
+++++++++++++++++++++++

By default (``neutralize: true``) the built system is made **net-neutral** with counter-ions, and ``salt_con`` can add extra salt.  How the ions are placed depends on the solvent:

- **Water** — VMD's ``autoionize`` replaces water molecules with ions (as usual).
- **Non-water** — ``autoionize`` cannot ionize a non-aqueous box (it needs water to make room), so pestifer ionizes by **replacing solvent molecules with ions**, a solvent-agnostic analogue of autoionize.

.. code-block:: yaml

   - solvate:
       solvent: DMSO
       # neutralize: true        # default; counter-ions cancel the solute's net charge
       # salt_con: 0.15          # (optional) also add this concentration of cation/anion pairs
       # cation: SOD             # defaults: cation SOD, anion CLA
       # anion: CLA

For the non-water path, ions are placed at the positions of solvent molecules chosen far from the solute and from one another (default 5 Å), which are deleted to make room, then an ``ION`` segment is built.  Ion valences are honored (e.g. a divalent ``CAL`` needs half as many to neutralize).  Set ``neutralize: false`` (and no ``salt_con``) to keep the system's net charge instead — NAMD then neutralizes it with a uniform background under PME.

