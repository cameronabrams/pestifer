.. _subs_buildtasks_ring_check:

ring_check
----------

The ``ring_check`` task detects and removes **pierced-ring** configurations — cases where a bond from one molecule passes through the geometric center of a ring in another molecule.  This most commonly arises in membrane systems containing cholesterol or other sterols, where a lipid acyl chain can thread through a sterol ring during initial packing.

Prerequisites
~~~~~~~~~~~~~

``ring_check`` requires PSF, PDB, and XSC state files, so it must follow at least one MD step (typically a short minimization).  The canonical placement is immediately after the first minimization following ``make_membrane_system``:

.. code-block:: yaml

   tasks:
     - make_membrane_system:
         ...
     - md:
         ensemble: minimize
         minimize: 1000
     - ring_check:

If no piercings are detected the task is a no-op, so it is safe to include unconditionally whenever sterols or glycans are present.

Parameters
~~~~~~~~~~

All parameters are optional:

.. code-block:: yaml

   - ring_check:
       segtypes: [lipid, glycan]   # segment types whose rings are checked (default: both)
       cutoff: 10.0                # bond–ring COM pre-screen distance in Å (default: 10.0)
       delete: piercee             # which molecule to remove (default: piercee)

**segtypes**
  List of segment types whose rings are examined.  Bonds from *any* segment type can be the piercer; only rings belonging to segments of these types are checked.  Default is ``[lipid, glycan]``.

**cutoff**
  Distance in Å between a bond midpoint and a ring center-of-mass used to pre-screen candidate pairs before the full geometric test.  Increasing it finds more candidates at the cost of more computation; the default of 10.0 Å is appropriate for typical ring sizes.

**delete**
  Controls which molecule is removed when a piercing is found:

  - ``piercee`` *(default)* — removes the ring-containing molecule (e.g. the cholesterol whose ring was threaded).  Use this for sterol rings pierced by lipid acyl chains.
  - ``piercer`` — removes the molecule whose bond did the piercing.  More appropriate when a glycan ring is pierced by a lipid chain and you want to keep the glycan.
  - ``none`` — logs all piercings but takes no action.  Useful for auditing.

How the detection works
~~~~~~~~~~~~~~~~~~~~~~~

For each ring in the target segments, pestifer uses a link-cell grid to find all bonds whose midpoint lies within ``cutoff`` Å of the ring's center of mass.  For each candidate bond, it projects the ring atoms onto a plane perpendicular to the bond at the bond midpoint and sums the subtended angles.  A sum near 2π indicates the bond passes through the ring interior (the winding-number test).

When piercings are found, the task writes a new psfgen script that deletes the offending residues and rebuilds the PSF/PDB, updating the pipeline state.

Example: sterol ring check in a viral membrane
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is taken from :ref:`example mper-tm viral bilayer` (Example 17):

.. code-block:: yaml

   tasks:
     - make_membrane_system:
         bilayer:
           composition:
             upper: {POPC: 0.5, POPE: 0.3, CHL1: 0.2}
             lower: {POPC: 0.5, POPE: 0.3, CHL1: 0.2}
     - md:
         ensemble: minimize
         minimize: 1000
         constraints:
           k: 10
           atoms: protein and name CA
     - ring_check:
