.. _subs_buildtasks_ligate:

ligate
------

A ``ligate`` task builds missing protein loops that are specified in the input structure file.  It is very often the case that a structure file declares "missing" or "zero-occupancy" residues that are part of the physical molecule that was analyzed but for whatever reason are not resolved.  Pestifer by default builds these in as unstructured loops and closes each internal gap.

By default (``method: ccd``) each internal loop is closed **offline and deterministically**:

1. The missing residues are inserted by the ``psfgen`` ``guesscoords`` command, which builds them from CHARMM internal coordinates as an essentially straight chain (every :math:`\phi`/:math:`\psi` = 180°) hanging off the last resolved residue, its far end *not* yet bonded to the downstream anchor.
2. Each loop residue's backbone dihedrals :math:`(\phi, \psi)` are re-seeded from the coil regions of the Ramachandran plot, giving a realistic, *compact* backbone rather than a straight chain.
3. The loop is closed onto the downstream anchor by **cyclic coordinate descent (CCD)**: the closure is distributed across the loop's own backbone torsions, rather than dragging one end across space.
4. An **iterative declash** refinement (perturb-and-reclose simulated annealing) removes steric clashes against the *actual* surroundings, including any loop copies already closed in this pass.  Each copy is closed **independently**, so copies converging at an assembly axis (e.g. the gp41 HR1N loops at a trimer 3-fold) interleave and avoid one another rather than being forced into an identical, mutually clashing conformation.
5. A final ``psfgen`` run applies the ``LINK`` patch to form the peptide bond between the loop and its downstream anchor (the ``connect`` step).

Because the closure is geometric, a few *soft* (non-threaded) steric overlaps may remain; the ``minimize`` task that ordinarily follows ``ligate`` relaxes these.  A loop whose backbone ends up *threaded* (a topological defect that minimization cannot undo) is flagged in the log, and the build can be made to abort on it via ``ccd.on_clash: error``.

Two optional knobs tune the closure: ``ccd.ensemble`` (number of independent Ramachandran-seeded starts per loop) and ``ccd.refine`` (iterations of iterative declash).  A ``ccd.guard: true`` option keeps the surrounding structure present *during* the CCD closure and rejects any rotation that increases the clash count, reaching fewer residual overlaps at the cost of a much slower closure; it is off by default because the downstream ``minimize`` does the same job.  ``ccd.thread_ca`` is the CA-displacement threshold (in Angstrom) that triggers a sequential re-close of any pair of loops whose backbones come closer than this after the parallel independent closures (default **0**, which is off).

The legacy ``method: steer`` performs a steered MD run that drags each loop's C-terminus onto its downstream anchor before the ``LINK`` patch is applied; it remains available as a fallback.  For insertions of two or fewer residues no closure is needed.

This is illustrated in :ref:`example env 4zmj`.
