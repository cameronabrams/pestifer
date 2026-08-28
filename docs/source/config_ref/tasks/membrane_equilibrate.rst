.. _config_ref tasks membrane_equilibrate:

``membrane_equilibrate``
========================

Self-terminating NPgT membrane equilibration -- the anisotropic sibling of density_equilibrate. Runs NPgT (useFlexibleCell, tensionless) in the same stability-bounded chunks and stops when BOTH the box density AND the membrane lateral area (A = a_x*b_y) have converged (density + area fully pin the anisotropic cell since V = A*c_z). By default a two-stage protocol (two_stage) first settles density at constant lateral area, then relaxes the area tensionless, so a freshly gridded (under-dense-in-z) bilayer does not co-compact its area down to a metastable gel-like APL. Convergence uses the same autocorrelation-corrected criterion per observable; area reuses the density tolerances unless area_* overrides are given. See docs/design/membrane-equilibrate.md.

Single-valued attributes:

  * ``temperature``: Temperature (K) (default: 300)

  * ``pressure``: Pressure (bar); NPT requires a solvated (periodic) system (default: 1)

  * ``two_stage``: If true (default), run the two-stage protocol: stage 1 settles the box density at CONSTANT lateral area (useConstantArea -- excess volume leaves from z only, so the build-correct lateral area is not dragged down) gating on density alone, then stage 2 relaxes the lateral area under tensionless NPgT (useConstantRatio) from the already-correct density, gating on both. Set false for the legacy single-stage tensionless run (density + area jointly), appropriate for a post-embed stage already near density equilibrium. (default: True)

  * ``constant_ratio``: If true (default), the tensionless area stage locks Lx:Ly (useConstantRatio) -- correct for an embedded membrane whose lateral aspect must stay stable. Set false to let the fully flexible cell relax x and y independently, appropriate for a calibration patch (no protein, no preferred aspect): on a non-square box (e.g. orthohexagonal) a locked skewed ratio makes the tensionless relaxation over-condense along the short axis instead of settling. (default: True)

  * ``dcdfreq``: number of time steps between configuration output to DCD file (default: 100)

  * ``xstfreq``: number of time steps between cell-size (extended-system) output; sets the density-sampling resolution (default: 100)

  * ``chunk_min``: nominal first-chunk length (steps) for the fast-shrinking start; crash recovery may temporarily go below this (default: 500)

  * ``chunk_max``: longest NPT chunk (steps); used once the box has settled (default: 5000)

  * ``shrink_safety``: size each chunk so the projected per-dimension cell shrink stays below this fraction of margin (default: 0.5)

  * ``margin``: NAMD patch margin (angstrom) -- the per-run cell-shrink slack used to bound chunk length (default: 4)

  * ``chunk_growth``: max factor by which a chunk may grow over the previous one (AIMD-style gentle re-growth after a crash-triggered halving) (default: 1.5)

  * ``max_shrink_retries``: consecutive patch-grid crash retries allowed at one boundary before giving up (each halves the chunk); the reactive stability net (default: 5)

  * ``max_steps``: hard ceiling on total NPgT steps; the task stops (with a warning) even if not converged. Higher than the isotropic density_equilibrate default because membrane post-embed densification (void removal) is inherently slow -- the hand-written ladder it replaces budgeted ~380k steps (default: 500000)

  * ``min_steps``: never declare convergence before this many total steps (default: 4000)

  * ``burn_in``: leading steps discarded (barostat transient) before assessing the density trend, measured relative to this task's first sample (not an absolute timestep) (default: 2000)

  * ``window_frac``: assess only the trailing this-fraction of post-burn-in samples, so the early densification ramp ages out of the window as the run lengthens (a full-history window never converges on a box that has plateaued) (default: 0.5)

  * ``drift_tol``: converged when the fractional DENSITY drift over the window falls below this (adequacy choice, ~0.2%; validated to converge robustly on a small/noisy box -- 0.1% is unreachable there in reasonable time) (default: 0.002)

  * ``precision_p``: precision gate -- require the (autocorrelation-corrected) SEM/mean of the windowed density below drift_tol/precision_p before testing drift (default: 3)

  * ``autocorr_reliability``: require the trailing window to span at least this many integrated autocorrelation times before its SEM is trusted (the density is autocorrelated over hundreds of steps, so a short window cannot resolve the mean no matter how flat it looks); the honest SEM = sigma/sqrt(N_eff), N_eff = N_window/tau_int, is size-aware for free since sigma/mean ~ 1/sqrt(N_atoms) (default: 6)

  * ``drift_conf``: sigmas added to the drift point estimate to form its upper confidence bound; convergence requires that bound to fall below drift_tol. Default 0 uses the point estimate -- the autocorrelation-corrected precision gate already prevents premature/spurious passes (a small point-drift cannot converge until the mean is resolved), so the bound is redundant and, because the trend is intrinsically ~sqrt(12) harder to bound than the mean, raising it sharply lengthens (or ceilings) small-box runs. Raise to ~1-2 only for an explicit statistical trend certificate. The slope SE is autocorrelation-inflated by sqrt(tau_int) (default: 0.0)

  * ``n_consecutive``: successive passing checks required before stopping (hysteresis against the noisy drift metric) (default: 3)

  * ``lipids_per_leaflet``: lipids per leaflet, used only to report area-per-lipid (APL = area / lipids_per_leaflet) in the report/plot; 0 = unknown (report raw area only) (default: 0)

  * ``area_drift_tol``: per-observable drift tolerance for the lateral AREA (soft/slow mode; defaults looser than density's, sized to the area's own finite-sample drift noise floor). Unset falls back to density drift_tol. (default: 0.012)

  * ``area_precision_p``: per-observable precision_p for the lateral area (relaxed vs density so the SEM gate = area_drift_tol/area_precision_p is attainable given the area's large autocorrelation time). Unset falls back to density. (default: 2.0)

  * ``area_drift_conf``: optional per-observable override of drift_conf for the lateral area (defaults to the density value)

  * ``area_burn_in``: per-observable burn_in for the lateral area, in steps (discards the initial fast area-relaxation transient before the window is evaluated). Unset falls back to density. (default: 20000)

  * ``area_min_steps``: minimum total steps before the AREA may be declared converged (a floor, ~1-2x the area autocorrelation time, that prevents premature convergence on a momentary flat spot during the initial relaxation). Sized for the 100/50 barostat (area tau ~2k steps); relax/tighten with the measured tau. Unset falls back to the shared min_steps. (default: 40000)

  * ``area_window_frac``: optional per-observable override of window_frac for the lateral area (defaults to the density value)

  * ``area_autocorr_reliability``: per-observable autocorr_reliability for the lateral area (relaxed vs density's 6 so the reliability guard is attainable given the long area autocorrelation time). Unset falls back to density. (default: 4.0)

  * ``area_plateau_tol``: optional cumulative area-plateau gate (fraction). When > 0, convergence also requires the lateral area to have genuinely flattened -- the quarter-over-quarter drift of the stage-2 area series (mean over the final quarter vs the preceding quarter) must be below this -- not merely the local windowed slope. This is the exact metric make_membrane_system uses to judge a calibration patch's preferred APL reliable (_AREA_CONVERGENCE_TOL), so set it (e.g. 0.005) on calibration patches to make "run until the creep flattens" the actual gate for a slowly-condensing (e.g. cholesterol-rich) leaflet. 0 (default) disables the extra gate. (default: 0.0)

  * ``other_parameters``: key:value pairs for other namd configuration file statements

  * ``single-core``: whether or not to run namd on a single core (default: False)

  * ``cpu-override``: revert to cpu version of namd3 (default: False)



Container-like attribute:

.. toctree::
   :maxdepth: 1

   membrane_equilibrate/addl_paramfiles


Subattribute:

.. toctree::
   :maxdepth: 1

   membrane_equilibrate/constraints


.. raw:: html

   <div class="autogen-footer">
     <p>This page was generated by ycleptic v2.3.0 on 2026-08-27.</p>
   </div>