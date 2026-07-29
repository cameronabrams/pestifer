# Design: `membrane_equilibrate` task (density + lateral-area convergence)

Status: **planned**. A self-terminating NPgT equilibration for membrane systems, the anisotropic
sibling of [`density_equilibrate`](density-equilibrate.md). It runs NPgT in stability-bounded chunks
and stops when **both** the box density **and** the membrane lateral area have converged, replacing the
hand-tuned fixed NPgT ladders used at three points in a membrane build. Reuses the
autocorrelation-corrected convergence machinery of `util/density_convergence.py` unchanged; the only
new physics is a second observable.

## Motivation — three fixed ladders, one engine

A membrane build hand-tunes a fixed equilibration ladder in three places, all of which want
self-termination (analogous to how `density_equilibrate` replaced the isotropic NPT ladder):

1. **Calibration-patch equilibration** (asymmetric builds, inside `make_membrane_system`) — relaxes
   small per-leaflet patches to measure preferred **area-per-lipid (APL)**, which *sizes the whole
   grid*. A poorly-converged APL mis-sizes the entire membrane, so rigorous area convergence matters
   most here. Today: a fixed NPgT stage + a crude two-block area-drift heuristic
   (`MakeMembraneSystemTask._area_convergence`, tol 0.5%).
2. **Bare-bilayer ("quilt") equilibration** (pre-embed, inside `make_membrane_system`) — fixed
   `minimize → NVT → NPT → NPgT 16k` ramp via `_autostage_protocol`.
3. **Post-embed equilibration** — a hand-written top-level `md:` ladder in examples 16/17:
   `NVT → NPT(2k/4k/8k/16k) → NPgT(16k/32k/64k/128k)`. The direct analog of the NPT ladder
   `density_equilibrate` replaced.

## Why density + area is the complete, minimal observable set

For an NPgT box with the membrane normal along **z**, the lateral area is `A = a_x·b_y` (or `|a×b|`
for a triclinic cell) and the volume is `V = A · c_z`. Hence

    density = mass / V = mass / (A · c_z).

If **A** and **density** both converge, then `c_z = mass/(density·A)` is determined — so the two
observables fully pin the anisotropic cell; thickness need not be tracked separately. Under NPgT with
`useConstantRatio` the x:y ratio is locked, so `A` is effectively a single soft lateral degree of
freedom relaxing under surface tension γ (default 0, tensionless).

## The engine generalizes almost for free

Everything hard in `util/density_convergence.py` — the autocorrelation-corrected SEM precision gate,
the drift test, the reliability guard, the shrink-rate chunk sizing, and the patch-grid crash/retry —
is already **observable-agnostic** (`DensityConvergenceMonitor` just consumes `(t, value)` samples).
Additions:

1. **`xst_cell_areas(path)`** — mirrors `xst_cell_volumes`; returns `(timesteps, areas)` with
   `A = a_x·b_y` (`|a×b|` for triclinic). `a_x`/`b_y` are already parsed by `_xst_data_rows`.
2. **`JointConvergence`** — a thin wrapper holding one monitor per named observable (density, area),
   feeding each its own series. A check passes only when **every** observable is stationary this
   chunk; a single shared `n_consecutive` counter means density and area must be *jointly* settled.
   Each observable keeps its own τ, SEM, and drift — they are independent statistically.
3. The chunk-stability loop is unchanged: `xst_max_shrink_rate` already keys off the max shrink across
   all three cell edges, so it correctly bounds an anisotropically condensing membrane; the reactive
   patch-grid crash/retry is orthogonal and reused verbatim.

## Task architecture — shared base + two thin subclasses

Refactor the shipped `density_equilibrate` into a small shared base carrying the chunk loop,
crash-retry, sizing, and reporting:

- `DensityEquilibrateTask` — `ensemble='NPT'`, observables `[density]` (unchanged, backward-compatible).
- `MembraneEquilibrateTask` — forces `ensemble='NPgT'` (so it picks up the `membrane` barostat block:
  `useFlexibleCell` + `useConstantRatio`, γ default 0), observables `[density, area]`.

Registered in `tasks/__init__.py`; contract `requires=(STATE,), provides=(STATE, MD_OUTPUT)` (same as
`MDTask`). The fuller "one `equilibrate` task with selectable `observables:`/`ensemble:`" (roadmap P3)
is a later consolidation this reaches incrementally.

## Tolerances — per-observable, defaulted equal, then **measured** for area

Area under γ=0 is a *soft* mode: its fluctuations are governed by the area-compressibility modulus
`K_A` and it plausibly has a **longer autocorrelation time and larger σ/mean than density** (membranes
relax and undulate slowly). So area may need its own tolerances — but this is an empirical question,
answered the same way we answered it for density, not a guess:

- Make `drift_tol`, `precision_p`, `drift_conf`, `autocorr_reliability`, and `burn_in`
  **per-observable** (area may override density's); default area's to density's values.
- **Measured (M2).** A fine-sampling (`xstfreq=1`) NPgT probe on the equilibrated 109k-atom DMPC
  membrane (example 16's final state) gives, on the plateau: **area** σ/mean ≈ 1.9e-3, **τ_int ≈ 921
  steps**; membrane **density** σ/mean ≈ 1.4e-3, τ_int ≈ 822 steps (cf. BPTI water density τ ≈ 559).
  So area *is* the slower, noisier mode — but only ~1.1–1.6× density's, the same order of magnitude.
  **Conclusion:** the density-defaulted tolerances are appropriate; area is simply the *binding*
  observable (it converges last), which is physically correct — the lateral area is the membrane's slow
  mode. No distinct area default is warranted; the `area_*` overrides remain available for a membrane
  that wants a deliberately looser area settling.

## Reporting / plotting

- `<basename>-density.dat` gains per-observable columns (density and area each with
  drift/drift_hi/SEM/τ/N_eff), plus **APL = area / lipids-per-leaflet** (lipid count from the PSF
  segtype or the `Bilayer` object).
- The plot becomes two-panel: density-vs-time and area-vs-time (APL on a twin axis), each with its own
  burn-in region, trailing window, and convergence marker.

## Integration — one engine, three sites

- **Post-embed (M3, easiest):** replace the hand-written NPgT ladders in examples 16/17 with a single
  `membrane_equilibrate:` task (keep the restrained minimize/NVT warm-up; replace the NPgT tail) — same
  migration pattern as the P2 protein work.
- **Bare-quilt + calibration-patch (M4, deeper):** `make_membrane_system` already runs its pre-embed
  relaxation through a subcontroller of nested tasks, so it can append a `membrane_equilibrate` task
  in place of `_autostage_protocol`'s fixed NPgT stages, and **retire the crude `_area_convergence`
  heuristic** for the rigorous monitor. Highest value on the calibration patch (it sizes the grid).

## Phasing (each shippable)

- **M1 — engine primitives (done).** `xst_cell_areas` (lateral area `|a×b|` from the `.xst`) and
  `JointConvergence` (converge only when all named observables are jointly stationary for
  `n_consecutive` checks; each sub-monitor keeps its own tolerances/τ) landed in
  `util/density_convergence.py`, NAMD-free and unit-tested. The **base-class refactor of
  `density_equilibrate` is deferred to M2**, where `MembraneEquilibrateTask` provides a concrete second
  consumer to design the shared base against — cleaner and lower-risk than a speculative refactor of the
  freshly-released task.
- **M2 — `membrane_equilibrate` task (done).** NPgT, density+area via `JointConvergence`, two-panel
  report/plot + APL (`lipids_per_leaflet`). The shared chunk loop was extracted into
  `ChunkedEquilibrateTask` (`equilibrate_base.py`); `density_equilibrate` and `membrane_equilibrate`
  are thin subclasses (density behavior verified unchanged by an end-to-end bpti1 build). Area
  tolerances default to the density values, with optional `area_*` per-observable overrides. Registered
  + schema-validated; 36 unit tests. **Area-autocorrelation calibration** (a fine-sampling `xstfreq=1`
  NPgT probe on a real DMPC membrane) measures whether area needs its own tolerances; end-to-end
  validation of the task in a full membrane build comes with M3.
- **M3 — migrate examples 16/17 post-embed ladders** to `membrane_equilibrate`.
- **M4 — wire the engine into `make_membrane_system` pre-embed** (quilt + calibration patch); retire
  the fixed NPgT tail of `_autostage_protocol` and `_area_convergence`.
- **Validation.** Full 16/17 builds; compare converged area/APL/density to the old fixed-ladder values;
  confirm APL calibration stability.

## Decisions (resolved)

- **Architecture:** shared-base refactor + `membrane_equilibrate` sibling (not the full unified
  `equilibrate` task yet).
- **Scope order:** M1–M3 (post-embed, lower-risk) first; M4 (`make_membrane_system` surgery) follows.
- **Surface tension:** γ=0 tensionless default (inherited from the `membrane` block); expose a `gamma`
  task spec as a documented hook but do not build the differential-stress workflow now.
- **APL:** report both raw area and APL.
- **Area tolerances:** per-observable, defaulted to density's, then set from the measured area
  autocorrelation in M2 (open until measured).
