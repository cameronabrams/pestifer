.. _subs_buildtasks_density_equilibrate:

density_equilibrate
-------------------

A ``density_equilibrate`` task is a **self-terminating** replacement for the hand-written ladder of
progressively longer NPT runs that used to end every solvated build.  It runs NPT and **stops itself**
when the box density has converged, or when a hard step ceiling (``max_steps``) is reached:

.. code-block:: yaml

   tasks:
     - ... (fetch / psfgen / solvate / minimize / NVT warm-up here)
     - density_equilibrate:
     - terminate:
         ...

With no attributes it uses the calibrated defaults; you rarely need to set anything.  It replaces a
block that previously looked like:

.. code-block:: yaml

   - md: {ensemble: NPT, nsteps: 200}
   - md: {ensemble: NPT, nsteps: 400}
   - md: {ensemble: NPT, nsteps: 800}
   - md: {ensemble: NPT, nsteps: 1600}
   - md: {ensemble: NPT, nsteps: 13200}
   - mdplot: {timeseries: [density], ...}

Why a *series of short* NPT runs
================================

NPT on a freshly solvated box shrinks the cell (the water was placed below its equilibrium density).
NAMD fixes its spatial decomposition (patch grid) and PME grid at the **start** of each ``run``; as the
cell shrinks *within* a run the patches shrink with it, and once a patch dimension falls below the
cutoff NAMD aborts (*"periodic cell has become too small for the patch grid"*).  Each new ``run`` (a
restart) recomputes the decomposition for the current, smaller cell.  So chunking is **mandatory for
numerical stability**, not merely a convergence convenience -- a single long NPT run on a shrinking box
*crashes*.

``density_equilibrate`` sizes each chunk from the shrink rate it just observed, so that the projected
per-dimension cell shrink stays below ``shrink_safety * margin``: chunks are short while the box shrinks
fast and lengthen as it settles -- the ladder's shape, derived rather than guessed.

How it decides it is done
=========================

Box density fluctuates intrinsically in the NPT ensemble (a small box legitimately fluctuates a few
tenths of a percent), and the fluctuation amplitude is size-dependent (``σ_ρ/mean ∝ 1/√N_atoms``).  The
task therefore tests the **trend against the statistical noise**:

#. discard the first ``burn_in`` steps (barostat transient);
#. a **precision gate** grows the averaging window (by running more chunks) until the standard error of
   the mean density is well below the tolerance (``SEM/mean < drift_tol / precision_p``) -- a small,
   noisy box automatically runs longer, a large quiet box stops sooner, so the *duration* (not the
   tolerance) carries the size dependence;
#. the box is **converged** when the fractional density drift over the window is below ``drift_tol`` for
   ``n_consecutive`` successive checks (hysteresis against the noisy drift metric).

``drift_tol`` (default ``1e-3``, i.e. ~0.1 % density stability) is a practical adequacy choice, not a
measured noise floor.  See ``docs/design/density-equilibrate.md`` for the full derivation and the
prototyping that rejected the alternatives (a fluctuation-amplitude threshold; a ``trend/SEM``
significance test).

Outputs
=======

Besides the usual continuation state (``psf``/``pdb``/``coor``/``vel``/``xsc``) that the next task
consumes, the task writes:

- a **convergence report** (``<basename>-density.dat``) -- per-chunk mean density, drift, ``SEM/mean``,
  precision-gate status, and the stop reason;
- a **density-vs-time plot** (``<basename>-density.png``).

If ``max_steps`` is reached without convergence the task stops anyway with a warning reporting the
residual drift and whether the precision gate was even met (a large residual, or an unmet gate, flags a
system that never settled).  A non-finite density (box blow-up) aborts immediately.

Scope
=====

This equilibrates the **box density** -- the job the NPT ladder did -- not full conformational
equilibration; a production run still follows.  It is for isotropic **solvated** boxes; membrane systems
should keep the ``NPgT`` protocol (anisotropic membrane equilibration is out of scope).

Because it wraps NAMD, the exact stop step is reproducible on a fixed machine/NAMD build with a fixed
seed, but -- unlike pestifer's offline numpy steps -- it is not bit-identical across platforms (GPU vs
CPU, NAMD versions); the *result* is an equilibrated box either way.

See the :ref:`config_ref tasks density_equilibrate` reference for the full attribute list.
