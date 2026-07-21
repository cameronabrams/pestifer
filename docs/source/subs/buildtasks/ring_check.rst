.. _subs_buildtasks_ring_check:

ring_check
----------

The ``ring_check`` task detects and resolves **pierced-ring** configurations — cases where a bond from one molecule passes through the geometric center of a ring in another molecule.  This most commonly arises in membrane systems containing cholesterol or other sterols, where a lipid acyl chain can thread through a sterol ring during initial packing.  It also happens with glycans, which contain multiple rings that can be threaded by (or thread through) protein side chains or other glycans.  Because these configurations are non-physical and destabilize downstream MD (they typically RATTLE-fail on the first dynamics step), pestifer provides this task to identify and resolve them.

How a piercing is resolved depends on the **segtype** of the pierced ring: protein and glycan rings are cleared by *rotation* (they cannot be deleted without breaking the molecule), while lipid rings are cleared by *deletion*.

Prerequisites
~~~~~~~~~~~~~

``ring_check`` requires PSF and PDB state files.  It is normally run after a minimization (which settles coordinates but cannot itself undo a topological piercing).  An XSC file is used when present: it supplies the periodic box vectors so the minimum-image convention can be applied when measuring bond–ring distances across periodic boundaries.  If no XSC is available (vacuum system), ``ring_check`` falls back to a non-periodic mode in which the bounding box is derived from the coordinate extents.

A common placement is immediately after the first minimization following a ``psfgen`` or ``make_membrane_system`` task.  (Note that a grid-packed ``make_membrane_system`` build now inserts its own lipid ``ring_check`` + minimize automatically before its first dynamics stage, so an explicit lipid ``ring_check`` there is mainly for belt-and-suspenders or for glycan/protein rings.)

Similarly, the ``psfgen`` task now resolves **glycan** ring piercings at graft time (rotating the offending glycan sub-branch out of the ring; see :ref:`source.sequence.glycans.declash.check_piercings <subs_buildtasks_psfgen>`), so an explicit glycan ``ring_check`` afterward is usually unnecessary — it remains useful as a safety check, or to resolve an aromatic protein ring by rotating its side chain (a motion the graft-time glycan-only rotation does not attempt) when the glycan rotation could not clear it.

.. code-block:: yaml

   tasks:
     - make_membrane_system:
         ...
     - md:
         ensemble: minimize
         minimize: 1000
     - ring_check:
         segtypes: [lipid, glycan]

Parameters
~~~~~~~~~~

.. code-block:: yaml

   - ring_check:
       segtypes: [lipid, glycan]   # REQUIRED: segment types whose rings are checked
       cutoff: 4.0                 # bond-ring COM pre-screen distance in Å (default: 4.0)
       max_ring_size: 7            # largest ring (atoms) to consider (default: 7)
       delete: piercee             # lipid-ring deletion strategy (default: piercee)

**segtypes**
  List of segment types whose rings are examined (e.g. ``lipid``, ``glycan``, ``protein``).  Bonds from *any* segment can be the piercer; only rings belonging to segments of these types are tested.  **There is no default** — if ``segtypes`` is omitted or empty, the task warns and checks nothing, so it must be set explicitly.

**cutoff**
  Distance in Å between a bond midpoint and a ring center-of-mass used to pre-screen candidate pairs before the full geometric (winding-number) test.  The default of 4.0 Å sits just above the 3.5 Å ring-piercing detection gate, so it prunes candidate bonds cheaply while still catching every true piercing.  A *larger* cutoff does not find more piercings; it only tests more candidate bonds and so runs slower on large membranes.

**max_ring_size**
  Only chordless cycles of at most this many atoms are treated as rings (default 7).  This bounds the cycle search so real 5-/6-membered chemical rings are found while the giant cycle that a disulfide closes through the protein backbone is skipped.

**delete**
  For a pierced **lipid** ring, which segment(s) to remove (this does not apply to protein or glycan rings, which are rotated):

  - ``piercee`` *(default)* — remove the ring-containing residue (e.g. the cholesterol whose ring was threaded).
  - ``piercer`` — remove the residue whose bond did the threading.
  - ``both`` — remove the piercee and, when it too is a lipid, the piercer.  A lipid tail threaded through a sterol ring strains *both* molecules, so removing both is the robust choice for inter-threaded lipids.
  - ``none`` — log all piercings and terminate the build.

.. note::

   **Resolution by segtype.**  A pierced **aromatic** protein side-chain ring (His/Phe/Tyr/Trp) is rotated off the piercing bond; a **rigid** ring that cannot itself move — a **proline** side-chain ring or a **glycan** ring — when speared by a glycan bond is cleared by rotating the glycan branch that carries the piercing bond.  A pierced **lipid** ring is deleted (per ``delete``).  If no available rotation clears a protein or glycan ring, the build stops with the offending residue named (rebuild, e.g. with a different random seed or more minimization).

How the detection works
~~~~~~~~~~~~~~~~~~~~~~~

For each ring in the target segments, pestifer uses a link-cell grid to find all bonds whose midpoint lies within ``cutoff`` Å of the ring's center of mass.  For each candidate bond, it projects the ring atoms onto a plane perpendicular to the bond at the bond midpoint and sums the subtended angles.  A sum near 2π indicates the bond passes through the ring interior (the winding-number test).

.. figure:: figs/ring_check_geometry.png
   :width: 80%
   :align: center
   :alt: Side view of a bond piercing a ring (left) vs. passing beside it (right).

   **Side view.**  Left: a bond (red) passes through the center of a hexagonal ring — a pierced configuration.  Right: the bond (blue) clears the ring.

.. figure:: figs/ring_check_winding.png
   :width: 90%
   :align: center
   :alt: Top-down view along the bond showing the winding-number test.

   **Winding-number test** (view along the bond axis).  Ring atoms are projected onto the plane perpendicular to the bond at the bond midpoint (star).  Left: the midpoint lies inside the ring; the angles subtended by consecutive atom pairs sum to ≈ 2π — pierced.  Right: the midpoint is outside the ring; the total subtended angle is much less than 2π — not pierced.

How pierced-ring resolution works
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A topological piercing cannot be undone by energy minimization: any continuous path that moves the bond out of the ring interior must cross the ring plane, and the force field has no term that penalizes the entanglement, so minimization settles the strained geometry without un-threading it.  ``ring_check`` therefore resolves piercings by discrete moves chosen per segtype:

- **Aromatic protein rings** (His/Phe/Tyr/Trp) pivot freely on their side-chain dihedrals, so the ring is swung off the piercing bond by rotating the side chain — ``chi2`` first, then ``chi1`` — trying a series of angles and re-checking after each.

- **Rigid rings** — a **proline** side-chain ring (fused to the backbone) or a **glycan** ring — cannot rotate themselves.  When such a ring is speared by a **glycan** bond, the glycan sub-branch carrying that bond is rotated as a rigid body about an upstream rotatable bond (a *bridge* of the glycan graph), pulling the bond out of the ring.

  Among the rotations that un-thread a ring, the one that introduces the fewest new close atomic contacts is kept, so a piercing fix never trades a thread for a steric clash.

- **Lipid rings** are resolved by deleting a segment per the ``delete`` parameter: a lipid molecule is an independent PSF segment, so removing it leaves the rest of the topology intact.  The deletion writes a psfgen script that issues ``delatom`` commands, rebuilds the PSF/PDB, and updates the pipeline state.

If every available rotation for a pierced protein or glycan ring fails, the build stops with the residue named and the manual fix suggested.

Example: sterol and glycan ring check in a viral membrane
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is taken from :ref:`example mper-tm viral bilayer` (Example 17), whose bilayer contains cholesterol (sterol rings) and whose embedded protein is glycosylated:

.. code-block:: yaml

   tasks:
     - make_membrane_system:
         ...
     - md:
         ensemble: minimize
         minimize: 1000
         constraints:
           k: 10
           atoms: protein and name CA
     - ring_check:
         segtypes: [lipid, glycan]
