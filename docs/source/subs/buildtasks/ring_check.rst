.. _subs_buildtasks_ring_check:

ring_check
----------

The ``ring_check`` task finds **pierced rings** (a covalent bond passing through the interior of
a ring) and resolves them by moving atoms. A pierced ring is topologically wrong chemistry that
no amount of minimization or dynamics undoes, so pestifer treats one as a build failure.

This page is the task reference. For how detection and resolution work, with pictures and a
worked example, see the guide :ref:`pierced_rings`.

Prerequisites
~~~~~~~~~~~~~

``ring_check`` needs a PSF and PDB in the pipeline state. If an XSC is present its cell is used,
so a bond and a ring on opposite sides of a periodic boundary are compared correctly. Without
one, the check runs in non-periodic mode.

Place a ``ring_check`` after any step that displaces atoms without dynamics (an embed, a
graft, a loop closure, initial packing), **before** the minimization that follows it, and again
after that minimization. Any move a resolving check makes must itself be minimized before
constrained dynamics; the pipeline validator refuses a pipeline that skips that step. The
pattern, and the reasons for it, are in :ref:`pierced_rings_where`.

You often do not need to add one yourself:

- ``psfgen`` resolves piercings from grafted glycans and modeled tails as it builds them.
- ``ligate`` resolves piercings from closing a loop, before patching it to its anchor.
- ``make_membrane_system`` wraps every minimization in its relaxation protocols with checks.
- ``terminate`` refuses to package any structure that still contains a piercing.

What those do not cover is the task list after an embed, which is why examples 16 and 17 carry
their own checks there.

Parameters
~~~~~~~~~~

.. code-block:: yaml

   - ring_check:
       segtypes: [lipid, glycan, protein]   # default; a list you write replaces it
       resolve: true                        # default; false = detect only
       delete: none                         # default; lipid rings only
       cutoff: 4.0                          # default, in Å
       max_ring_size: 7                     # default

**segtypes**
  Segment types whose rings are examined. Bonds from any segment can be the piercer. The default
  checks lipid, glycan and protein rings. A list you write *replaces* the default rather than
  adding to it, so ``segtypes: [lipid]`` checks only lipid rings. An empty list is an error,
  since a check of nothing would report success.

**resolve**
  ``true`` (the default) resolves what it finds by moving atoms. ``false`` only detects, and any
  piercing found stops the build. Use ``false`` for the last check before dynamics, where a move
  would leave an unrelaxed structure.

**delete**
  What to do with a piercing of a **lipid** ring that no move clears: ``none`` (the default)
  stops the build and names every remaining piercing. ``piercee`` deletes the residue owning
  the ring, ``piercer`` the residue owning the bond, and ``both`` deletes both when both are
  lipids. Deleting changes the composition you asked for, so only use it when losing a lipid is
  preferable to losing the build. Protein and glycan rings are never deleted. A remaining
  piercing involving either is always fatal.

**cutoff**
  Radius in Å within which a bond's midpoint must lie from a ring's center to be tested. Any
  value at or above the piercing test's own reach gives identical results. Larger values only
  test more bonds and run slower on big membranes.

**max_ring_size**
  Largest ring, in atoms, that is considered. The default of 7 covers chemical rings and skips
  the large cycles a disulfide bond closes through a protein backbone.

Output
~~~~~~

The task records what it found and did: ``piercings_initial``, ``rotations`` and
``piercings_remaining`` appear in ``run-record.json``, so whether a pass moved anything is a
recorded fact rather than something to read out of the log. Every failure message contains the
literal token ``RING_PIERCING_FOUND``. Search for that token, not the surrounding wording.

Example: after embedding a protein in a sterol-rich bilayer
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From :ref:`example mper-tm viral bilayer`. The embed is checked on both sides of a minimization,
the second check's moves are relaxed by another minimization, and the last check only looks:

.. code-block:: yaml

   tasks:
     - make_membrane_system:
         ...
     - ring_check:
         segtypes: [lipid, glycan, protein]
     - md:
         ensemble: minimize
         minimize: 1000
         constraints:
           k: 10
           atoms: protein and name CA
     - ring_check:
         segtypes: [lipid, glycan, protein]
     - md:
         ensemble: minimize
         minimize: 1000
         constraints:
           k: 10
           atoms: protein and name CA
     - ring_check:
         segtypes: [lipid, glycan, protein]
         resolve: false
