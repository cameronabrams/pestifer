# Design: `density_equilibrate` task

Status: **planned** — this doc is the plan; nothing built yet. Direction: a single post-solvation
task that runs NPT until the box density has converged (or a step ceiling is hit), replacing the
fixed ladder of progressively longer NPT runs. Fully offline, reuses pestifer's existing MD /
continuation substrate.

## Problem

Every solvated build today ends with a hand-written ladder of NPT runs of increasing length, e.g.

```yaml
- md: {ensemble: NPT, nsteps: 200}
- md: {ensemble: NPT, nsteps: 400}
- md: {ensemble: NPT, nsteps: 800}
- md: {ensemble: NPT, nsteps: 1600}
- md: {ensemble: NPT, nsteps: 13200}
```

This is a **heuristic with magic numbers baked into every config**. It does not adapt to the
system: a small soluble peptide box over-equilibrates on this schedule while a large glycoprotein
trimer or a membrane patch under-equilibrates — both run the same fixed number of steps. It is also
noise: the same five lines are copied into the template, the `new-system` pipeline, and every
bundled example.

## Goal / non-goals

- **Goal:** one `density_equilibrate` task that runs NPT and **stops itself** when the system's box
  density has converged, or when a hard step ceiling is reached — replacing the ladder across the
  board (template, `new-system` interactive pipeline, bundled examples).
- **Non-goals (v1):**
  - *Full-system* (conformational) equilibration — density is a **box/solvent** observable (see
    below); protein RMSD / potential-energy relaxation is out of scope.
  - A general pluggable-observable equilibrator — v1 is density-only. Generalizing to an
    `equilibrate` task with a selectable observable is a possible P3.

## What "density converged" actually means

The NPT ladder exists to let the **solvent box relax to its equilibrium density** after the abrupt
insertion of a solute into a freshly built water box: the barostat contracts/expands the cell until
the water reaches ~1 g/cc at the target pressure. Density is therefore *exactly* the right
observable for the job this task replaces.

It is **not** a claim that the whole system is equilibrated. Box density plateaus much faster than
conformational relaxation (backbone RMSD, side-chain repacking, slow H-bond networks). The task is
named and documented for what it does — *equilibrate the box density* — not "the system is
equilibrated." A production run still follows.

## Approach — chunked NPT + convergence test

Run NPT in **chunks** of `chunk_steps`; after each chunk, extract the density time series so far and
test for convergence; stop when converged or `max_steps` is reached.

```
firsttimestep = 0
series = []                              # (timestep, density) samples
while total_steps < max_steps:
    run NPT for chunk_steps (continuation from the running state; firsttimestep = total_steps)
    series += density_samples_from_xst(this chunk)
    total_steps += chunk_steps
    if total_steps >= min_steps and converged(series):
        break
register the final state; write a convergence report + density-vs-time plot
```

This reuses the existing substrate directly:

- **Continuation.** `MDTask.namdrun` already threads `firsttimestep` from the state artifact and
  each run emits the state fileset (`psf`/`pdb`/`coor`/`xsc`/`vel`) that the next chunk continues
  from (`tasks/mdtask.py`, `tasks/continuation.py`). The chunk loop is just repeated `namdrun`
  calls with the running `firsttimestep`. This also dovetails with the roadmap's
  *restart/resume an interrupted build* item — same per-task state fileset.
- **Density.** NAMD writes the extended-system trajectory `.xst`, one line per `xstFreq` frame:
  `ts a_x a_y a_z  b_x b_y b_z  c_x c_y c_z  ...`. For an orthorhombic cell the volume is
  `a_x·b_y·c_z` (the same columns `_parse_xsc_cell` already reads in
  `util/densityprofile.py`: indices 1, 5, 9). Total system mass is known from the PSF. So
  `density(t) = M_total / (a_x·b_y·c_z)` per frame — a scalar density time series with no new
  machinery. (For a non-orthorhombic cell, volume is the scalar triple product of the three cell
  vectors.)

## The convergence criterion (the crux)

**The observable is equilibrated when its *mean* has stopped drifting — not when its fluctuations
are small.** Density fluctuations are *intrinsic to the NPT ensemble* and do **not** shrink as the
system equilibrates (a small box legitimately fluctuates ±0.3–0.5 %). So the criterion tests the
**trend against the statistical noise**, i.e. whether the systematic drift across the window is
smaller than the uncertainty in the mean.

*(This corrects an earlier version of this doc that also required the fluctuation amplitude itself
to fall below a threshold — validated wrong against real `.xst` runs: an equilibrated soluble box
never satisfies it, because its fluctuations are intrinsic. See "Empirical calibration" below.)*

**Fractional-plateau test with hysteresis.** (A `trend < k·SEM` *significance* test was prototyped
and **rejected** — see "Empirical calibration": the trend/SEM ratio is noisy, sampling-dependent,
and reads ~5–6 for every equilibrated run, so no universal `k` exists. Use a **scale-independent
fractional tolerance**, which clusters tightly across box sizes.)

1. **Burn-in.** Discard the first `burn_in` steps (the initial barostat transient).
2. **Window.** Consider the trailing `window` steps (or all post-burn-in samples if shorter).
3. **Block-average** into `n_blocks` contiguous blocks → block means `d̄_k` at `t_k`.
4. **Fractional drift.** `drift = |slope| · (t_last − t_first) / mean` from a line fit to
   `(t_k, d̄_k)` — the fractional change of the mean density across the window. (Equivalently, a
   two-window Δ: `|mean(2nd half) − mean(1st half)| / mean`; both track together.)
5. **Converged** when `drift < drift_tol` for **`n_consecutive` successive checks** (hysteresis),
   provided `total_steps ≥ min_steps`.

The **hysteresis is essential**: the drift metric is noisy and non-monotonic — even a still-creeping
box shows the occasional low-drift window (2bgj dips to `6.6e-4` at one check, then back to
`4.5e-3`), so a single passing check false-triggers. Requiring several consecutive passes stops on a
sustained plateau, not a lucky window.

If `max_steps` is reached without convergence, the task **stops anyway** with a loud warning
reporting the residual `drift` (a large residual flags a system that never settled), rather than
running forever. A NaN / box-blowup in the NAMD output aborts immediately.

**Tolerance is not yet pinned — it needs a proper reference run (open).** The archived `.xst` files
are all too short to establish it: they are 16k-step equilibrations whose density is *still creeping
~0.3–0.4 %* at their end (2bgj rises 1.002 → 1.005 from step 10k to 16k), or drifting membranes.
Before fixing `drift_tol`, generate one **long** NPT reference run on a small box (e.g. 100–200k
steps) so the density visibly plateaus, and read off the residual drift *at plateau* — that sets the
tolerance floor. Until then `drift_tol` is a placeholder.

## Empirical calibration (prototype on real `.xst` runs)

The criterion was prototyped against archived pestifer NPT `.xst` files (density taken as
`1/volume`, since the drift/noise metrics are all *fractional* — no PSF mass needed):

- **Soluble target (2bgj, a small solvated protein).** Density jumps 0.906 → ~1.00 within ~1400
  steps of the barostat contracting the box, then just **fluctuates** around ~1.005 for the
  remaining ~17 000 steps of the old ladder — i.e. it is equilibrated by ~step 3 000, and the ladder
  wastes the rest. This is the premise, confirmed. The trend-vs-noise test stops it at **~7–9k
  steps** — versus the ladder's ~20k — a real ~60 % saving, with the mean well-established.
- **Failure modes it exposed, in order.** (a) `fluct_tol` never lets 2bgj converge (intrinsic
  fluctuations) — removed. (b) A `trend < k·SEM` significance test looked promising on 2bgj alone
  (converged at `~1.5·SEM`, ~7–9k steps) but **does not generalize**: across five systems (2bgj,
  GroEL/1aon, 6cm3, an AlphaFold model, ex20) the trend/SEM ratio at *end-of-run* is a consistent
  **~5–6**, not ~2 — 2bgj's earlier "convergence" was a noise dip (back to 5.68 by step 16k). SEM is
  sampling-dependent, so no universal `k` exists — rejected. (c) The **scale-independent fractional**
  metrics, by contrast, cluster tightly across those same very-different box sizes:
  `trend/mean ≈ 2–4.5e-3`, two-window `Δ ≈ 1–3.5e-3`. That is the criterion to use.
- **The metric is noisy → hysteresis needed.** Tracing 2bgj's trailing-window `drift` over the run:
  it is 4.5e-2 → 2.7e-2 (equilibrating) then jumps 6.6e-4 → 3.8e-3 → 4.5e-3 around the plateau. A
  single-window threshold fires on the 6.6e-4 dip prematurely; `n_consecutive` passes fix this.
- **The archived runs are too short to set the tolerance.** They are 16k-step equilibrations still
  creeping ~0.3–0.4 % at their end (2bgj mean density 1.002 → 1.005 from step 10k → 16k) — not a true
  plateau. A dedicated long reference run is required before `drift_tol` is fixed (see criterion).
- **Scope note.** The membrane npat run drifts ~6 % monotonically over 256k steps and never
  plateaus — anisotropic membrane equilibration is **out of scope** (keep the npat/npgt protocol);
  the `max_steps` ceiling + residual-drift report are the backstop for a system that never settles.

Provisional defaults (form settled, tolerance pending calibration): `window ≈ 10000`,
`n_blocks ≈ 6`, `burn_in ≈ 2000`, `min_steps ≈ 4000`, `chunk_steps ≈ 2000`, `n_consecutive ≈ 3`;
`drift_tol` **TBD** from a long reference run (provisionally ~1–2e-3, i.e. hold the mean density
to ≲0.1–0.2 % across the window for several checks).

## Parameters (task spec)

| key            | meaning                                                          | default (calibrated) |
|----------------|------------------------------------------------------------------|----------------------|
| `chunk_steps`  | NPT steps per convergence check                                  | 2000                 |
| `max_steps`    | hard ceiling; stop even if not converged                         | 100000               |
| `min_steps`    | never declare convergence before this many steps                 | 4000                 |
| `burn_in`      | leading steps discarded before assessing the trend               | 2000                 |
| `window`       | trailing steps the trend/noise test is computed over             | 10000                |
| `n_blocks`     | blocks the window is averaged into                               | 6                    |
| `n_consecutive`| successive checks that must pass before stopping (hysteresis)    | 3                    |
| `drift_tol`    | converge when the fractional drift `trend/mean < drift_tol`      | **TBD** (~1–2e-3)    |
| `seed`         | NAMD RNG seed (reproducible stop on a fixed machine)             | fixed default        |

Plus the usual NPT knobs it forwards to `namdrun` (`temperature`, `pressure`, `timestep`, output
frequencies). Defaults are starting points to be tuned during validation (see below).

## Task integration

- New task, YAML header `density_equilibrate`, subclassing `MDTask` (it *is* an MD run, with a
  convergence loop around it). `pipeline_contract`: `requires=(STATE,)`, `provides=(STATE,)`.
- `do()` implements the chunk loop above, calling the existing `namdrun(ensemble='NPT', ...)` per
  chunk with the running `firsttimestep`, accumulating density from each chunk's `.xst`.
- Registers the final `state` fileset (as any MD task does) and two artifacts: a **convergence
  report** (per-chunk drift/fluct, the stop reason) and a **density-vs-time plot** (reuse the
  `mdplot` density machinery).
- Schema: a `density_equilibrate` task entry mirroring `md`'s NPT options plus the table above.

## Reproducibility

Pestifer values reproducible builds. With an explicit `seed`, NPT is deterministic on a **fixed
machine/NAMD build**, so the adaptive stop point is itself reproducible — the stop step is part of
the (seeded) trajectory. It is **not** bit-identical across platforms (GPU vs CPU, NAMD versions):
the trajectory, and hence the exact stop step, may differ. The *result* is an equilibrated box
either way, but — unlike the offline numpy work (CCD, declash) — this task is not cross-platform
bit-reproducible, by nature of wrapping NAMD. Documented as such.

## Validation (first-class deliverable)

**Prerequisite — calibrate the tolerance.** First run *one* long NPT reference (100–200k steps) on a
small solvated box until the density visibly plateaus; measure the residual `drift`/two-window Δ at
plateau; set `drift_tol` just above that floor. This is a hard prerequisite — the criterion form is
settled but the tolerance is not, and it cannot be read off the (too-short) archived runs.

Then, run on systems spanning the size range and compare against a long fixed NPT run:

- a small soluble protein (fast box, should stop early),
- a large glycoprotein trimer (e.g. an HIV-Env example),
- a membrane patch (anisotropic box; confirm the orthorhombic-volume density is sensible, or use
  the full triple-product volume).

Confirm: (a) the final density matches a long reference run within tolerance, (b) the stop step is
sensible (not premature, not the ceiling), (c) determinism on a fixed machine (same seed → same
stop). Tune the default tolerances/window from these runs.

## Phasing (each shippable)

- **P1 — the task.** `density_equilibrate` with the drift+fluctuation criterion, density from the
  `.xst`, the `max_steps` guard and NaN abort. Wire it into the template and the `new-system`
  interactive pipeline (replacing the NPT ladder there). Validate on the three systems above.
- **P2 — migrate the bundled examples** off the NPT ladder onto `density_equilibrate`; add the
  convergence-report + density-plot artifacts. Keep plain `md` available for users who want a fixed
  schedule.
- **P3 (optional) — generalize** to an `equilibrate` task with a selectable observable (density,
  box dimensions, potential energy, an arbitrary `.xst`/log column), reusing the same
  chunk-and-test loop.

## Decisions

- **Fixed:** chunked NPT with convergence-or-ceiling stopping; density from the `.xst` cell volume
  (triple product) + PSF mass; convergence by a **scale-independent fractional-drift tolerance with
  `n_consecutive`-check hysteresis** — *not* a fluctuation threshold and *not* a `trend/SEM`
  significance test (both prototyped and rejected against real runs); density-only for v1; a distinct
  task (not an `md` `until:` mode); seeded / deterministic-per-machine; isotropic **soluble-box**
  scope (membranes keep npat/npgt); replace the ladder across template + `new-system` + examples;
  keep plain `md` for fixed schedules.
- **Open:** `drift_tol` itself — **must be calibrated from a long reference run** (prerequisite in
  Validation), not from the too-short archived runs; final `window`/`n_consecutive`; the
  convergence-report/plot artifact format.
