.. _subs_buildtasks_merge:

merge
-----

A ``merge`` task combines two or more independently prepared PSF/PDB systems into a single, unified PSF/PDB system.  It is useful when you want to assemble a multi-component system (e.g. two protein copies, a protein and a nucleic acid, or any combination of separately built subsystems) without rebuilding the topology from scratch.

The task uses VMD and psfgen internally: it first renames any conflicting segment names across the input systems, then merges all systems into one via psfgen ``readpsf``.  Patch remarks (such as disulfide bond ``DISU`` patches) are preserved for all input systems, with segment names updated to reflect any renames.

Each entry in the ``systems`` list must provide a ``psf`` file and either a ``pdb`` or ``coor`` coordinate file.  An optional ``segname_map`` dict can supply explicit per-segment renames for that system.

Segment name collisions between systems are resolved automatically according to ``collision_strategy``:

* ``enumerate`` *(default)*: the conflicting name is shortened to three characters and a numeric suffix (``0``–``9``, then ``00``–``99``) is appended until a unique name is found.  For example, ``PROA`` becomes ``PRO0``, ``PRO1``, and so on.
* ``error``: a :class:`ValueError` is raised immediately on the first unresolved collision.  Use this when every rename must be explicit.

Example
+++++++

The following example merges two copies of a pre-built BPTI vacuum system.  The second copy is a PDB with all atom coordinates translated 100 Å along Z so the two molecules do not overlap.  Segment name collisions are resolved automatically by enumeration.

.. code-block:: yaml

   tasks:
     - merge:
         systems:
           - psf: bpti.psf
             pdb: bpti.pdb
           - psf: bpti.psf
             pdb: bpti-shifted.pdb
         collision_strategy: enumerate

After the task, the merged state (a single PSF and PDB) is available to downstream tasks such as :ref:`solvate <subs_buildtasks_solvate>` or :ref:`md <subs_buildtasks_mdtask>`.

Explicit segment renaming
+++++++++++++++++++++++++

If you prefer to control the final segment names yourself, supply a ``segname_map`` for the system whose segments need renaming:

.. code-block:: yaml

   tasks:
     - merge:
         systems:
           - psf: bpti.psf
             pdb: bpti.pdb
           - psf: bpti.psf
             pdb: bpti-shifted.pdb
             segname_map:
               A: AA
               B: BB
               C: CC

Explicit maps take priority over the automatic enumeration strategy.  If an explicitly requested name still collides with an earlier system, the automatic strategy is used as a fallback.

Checking the interfaces the merge creates
+++++++++++++++++++++++++++++++++++++++++

Each input was built and checked on its own, so the contacts *between* inputs are first seen here.  A part modeled in one input -- a loop the crystal structure did not resolve, say -- knows nothing about where another input sits, and can arrive inside it.  After merging, the task therefore checks every pair of atoms and every ring across the inputs:

* Heavy atoms from different inputs closer than 1.2 Å are reported at WARNING, by residue pair.  A few are what the minimize that should follow a merge is for; many, concentrated on one region, mean that region was built into space another input occupies.
* A ring (protein, glycan or lipid) pierced by a bond from a *different* input stops the build.  The inputs are relaxed structures, and a minimize started from a relaxed structure locks a threading in rather than undoing it, so ``terminate`` would refuse it anyway -- after every stage in between had run.

``interface_check`` sets the response to a piercing: ``error`` *(default)*, ``warn`` (report it and continue), or ``off`` (skip both checks).  The fix is to build the input that models the region with the other input present as an obstacle; see ``obstacles`` under :ref:`psfgen <subs_buildtasks_psfgen>`.

For a full description of all parameters see :ref:`config_ref tasks merge`.
