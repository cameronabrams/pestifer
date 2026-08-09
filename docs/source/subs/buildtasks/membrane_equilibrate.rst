.. _subs_buildtasks_membrane_equilibrate:

membrane_equilibrate
--------------------

A ``membrane_equilibrate`` task is the anisotropic sibling of
:ref:`density_equilibrate <subs_buildtasks_density_equilibrate>`: a **self-terminating** NPgT
equilibration that replaces the hand-tuned NPgT ladder a membrane build used to need.  It runs
tensionless NPgT (``useFlexibleCell``) in stability-bounded chunks and stops when **both** the box
density and the membrane's lateral area have converged:

.. code-block:: yaml

   tasks:
     - ... (make_membrane_system, which embeds the protein, here)
     - membrane_equilibrate:
     - terminate:
         ...

With no attributes it uses the calibrated defaults.  It replaces a post-embed block that previously
looked like:

.. code-block:: yaml

   - md: {ensemble: NVT,  nsteps: 1000}
   - md: {ensemble: NPT,  nsteps: 2000}
   - md: {ensemble: NPT,  nsteps: 4000}
   - md: {ensemble: NPT,  nsteps: 8000}
   - md: {ensemble: NPT,  nsteps: 16000}
   - md: {ensemble: NPgT, nsteps: 16000}
   - md: {ensemble: NPgT, nsteps: 32000}
   - md: {ensemble: NPgT, nsteps: 64000}
   - md: {ensemble: NPgT, nsteps: 128000}

The task also runs *inside* ``make_membrane_system``, as a stage of a ``relaxation_protocols``
protocol, where it relaxes the calibration patches and the bare quilt.  See
:ref:`make_membrane_system <subs_buildtasks_make_membrane_system>` for that usage; everything on this
page applies there too.

Why density **and** area
========================

With the membrane normal along :math:`z`, the lateral area is :math:`A = a_x b_y` and the cell volume
is :math:`V = A c_z`, so

.. math::

   \rho = \frac{m}{V} = \frac{m}{A\,c_z}.

If density and area have both converged, then :math:`c_z = m/(\rho A)` is determined as well — the two
observables fully pin the anisotropic cell, and the box thickness needs no separate criterion.  That is
what makes this the *minimal complete* observable set for a membrane, where tracking density alone
(as ``density_equilibrate`` does) would leave the lateral degree of freedom unconstrained.

Each observable is tested independently by the same autocorrelation-corrected machinery
``density_equilibrate`` uses — its own integrated autocorrelation time, its own SEM precision gate, its
own drift test — and a check passes only when *every* observable is stationary.  A single shared
``n_consecutive`` counter (default 3) means density and area must be **jointly** settled, not merely
each settled at some point.

The two-stage protocol
======================

A freshly gridded bilayer is built at roughly the right lateral area but is **under-dense in** :math:`z`,
because the per-leaflet z-reservation over-sizes the box height.  A single tensionless piston asked to
fix the density *and* find the equilibrium area at once takes the excess volume out of the lateral
dimensions too, co-compacting the area down to a metastable, gel-like area per lipid that subsequent MD
will not undo.

``membrane_equilibrate`` therefore decouples the two by default (``two_stage: true``):

**Stage 1 — settle the density at constant lateral area.**  ``useConstantArea`` freezes :math:`x` and
:math:`y`, so the excess volume can only leave from :math:`z` and the build-correct area is not dragged
down.  Only density is gated here; the area is pinned, so it is not yet a meaningful observable (the
report still prints the measured area and APL, with drift shown as ``n/a``).

**Stage 2 — relax the area, tensionless.**  Once density has settled, the barostat switches to
tensionless NPgT and the lateral area relaxes from an already-correct density — a small move that no
longer collapses.  Both observables are gated.  Both monitors are reset at the hand-off so stage-2
convergence is judged only on stage-2 samples.

Set ``two_stage: false`` for the legacy single-stage behavior, which is reasonable for a stage that is
already near density equilibrium.

.. note::

   ``constant_ratio`` (default ``true``) locks the :math:`L_x{:}L_y` ratio during the tensionless
   stage.  That is correct for an **embedded** membrane, whose lateral aspect must stay stable.  It is
   wrong for a **calibration patch**, which has no protein and no preferred aspect: on a non-square
   (orthohexagonal) box a locked skewed ratio makes the relaxation over-condense along the short axis
   instead of settling.  Set ``constant_ratio: false`` on calibration protocols.

Why the chunks
==============

For the same reason ``density_equilibrate`` chunks its NPT: NAMD fixes its spatial decomposition and PME
grid at the **start** of each ``run``, and if a cell dimension shrinks past the patch margin *within* a
run, NAMD aborts with *"periodic cell has become too small for the patch grid"*.  Each chunk is sized
from the shrink rate just observed so the projected shrink stays below ``shrink_safety * margin``, and a
patch-grid crash is caught and retried with a halved chunk (up to ``max_shrink_retries``).  The chunk
sizing keys off the largest shrink across all three cell edges, so it bounds an anisotropically
condensing membrane correctly.

Area per lipid
==============

Where it can, the task measures the **leaflet geometry** from the incoming frame — the lipid count in
each leaflet and the protein's cross-sectional footprint in each — and reports APL per leaflet and
**protein-corrected**:

.. math::

   \mathrm{APL} = \frac{A - A_\mathrm{protein}}{n_\mathrm{lipids}}

which is the number worth comparing against a literature APL.  If the geometry cannot be measured the
task warns and falls back to the naive ``area / lipids_per_leaflet`` using the ``lipids_per_leaflet``
spec, if you supplied one; with neither, it reports raw area only.

Tolerances
==========

Area is a **soft, slow** mode — it relaxes under the area-compressibility modulus and has a much longer
autocorrelation time than the density — so it carries its own, looser defaults rather than reusing the
density's: ``area_drift_tol`` (1.2 %, vs. 0.2 % for density), ``area_precision_p``,
``area_autocorr_reliability``, ``area_burn_in`` (20 000 steps, discarding the initial fast area
transient), and ``area_min_steps`` (40 000 steps, a floor of roughly one to two area autocorrelation
times that prevents certifying a momentary flat spot mid-relaxation).  Any ``area_*`` spec left unset
falls back to the corresponding density value.

``area_plateau_tol`` adds an optional, stricter gate on top of the local one.  The windowed-slope test is
*local*: a slowly-creeping cholesterol-rich leaflet can show a small instantaneous slope while a
systematic descent continues.  When ``area_plateau_tol`` is greater than zero, convergence additionally
requires the **cumulative** quarter-over-quarter drift of the stage-2 area series (mean of the final
quarter vs. the preceding quarter) to fall below it.  This is exactly the metric
``make_membrane_system`` uses to judge a calibration patch's APL reliable, so setting it there (e.g.
``0.005``) makes "run until the creep flattens" the actual stopping rule, and the downstream
reliability check then agrees with the task's own verdict.

Outputs
=======

Besides the continuation state (``psf``/``pdb``/``coor``/``vel``/``xsc``) the next task consumes, the
task writes:

- a **convergence report** (``<basename>-membrane.dat``) — per-chunk stage, density, area, per-leaflet
  APL, signed drifts, ``SEM/mean``, pass count and reason, with the leaflet geometry and the stage-1 →
  stage-2 hand-off step in the header;
- a **two-panel plot** (``mdplots/<basename>-membrane.png``) — density and lateral area vs. time, with a
  secondary APL axis, the burn-in and averaging windows shaded, and the stage hand-off and convergence
  steps marked.

Drifts are reported **signed** (negative = shrinking) so the direction of a trend is legible; the
convergence gate keys off the magnitude.

If ``max_steps`` is reached without convergence the task stops with a warning that reports the residual
drift of each observable and whether its precision gate was even met — a large residual, or an unmet
gate, flags a system that never settled.  Its default (500 000 steps) is higher than
``density_equilibrate``'s because post-embed densification, which has to squeeze out the voids left
around the protein, is inherently slow.  A non-finite density or area (box blow-up) aborts immediately.

See ``docs/design/membrane-equilibrate.md`` for the derivation and the
:ref:`config_ref tasks membrane_equilibrate` reference for the full attribute list.
