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

You must first **build and install** the solvent box (see :ref:`sub_make_pdb_collection`):

.. code-block:: bash

   $ pestifer make-pdb-collection solvent --resname MEOH --nmol 216 --density 0.792
   $ pestifer modify-package pdb-repo add-entry solvent/MEOH --collection solvent

.. note::

   Building a solvent box currently requires a coordinate source for the single molecule -- either an existing PDB-repository entry (as water and ions have) or a complete internal-coordinate (IC) table in the residue's CHARMM topology.  Many small-molecule CGenFF solvents (e.g. ``MEOH``, ``ETOH``, ``DMSO``) ship without an IC table, so a starting single-molecule PDB for them is not yet auto-generated.

