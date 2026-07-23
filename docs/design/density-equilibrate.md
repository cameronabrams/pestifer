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

**Trend-vs-noise test.**

1. **Burn-in.** Discard the first `burn_in` steps (the initial barostat transient).
2. **Window.** Consider the trailing `window` steps (or all post-burn-in samples if shorter).
3. **Block-average.** Split the window into `n_blocks` contiguous blocks; take each block's mean
   density `d̄_k` at its mid-time `t_k`. Blocks larger than the density autocorrelation time make
   the `d̄_k` ~independent, so their spread estimates the sampling error of the mean.
4. **Trend.** Least-squares-fit a line to `(t_k, d̄_k)`; the systematic change across the window is
   `trend = |slope| · (t_last − t_first)`.
5. **Noise.** The standard error of the mean is `SEM = std(d̄_k) / √n_blocks`.
6. **Converged** when `trend < k · SEM` (the trend is not significant against the sampling noise),
   provided `total_steps ≥ min_steps`. A secondary absolute ceiling `trend / mean < drift_tol`
   guards the degenerate case of a window too noisy to resolve any trend.

Testing the trend *relative to the noise* is what lets an equilibrated-but-fluctuating box converge
(2bgj: `trend ≈ 1.5·SEM` at ~7–9k steps) while a genuinely creeping system does not (a membrane
npat run drifting ~6 % monotonically stays `trend ≫ SEM` far longer).

If `max_steps` is reached without convergence, the task **stops anyway** with a loud warning that
reports the residual `trend`/`SEM` (a large residual flags a system that never settled — worth the
user's attention), rather than running forever. A NaN / box-blowup in the NAMD output aborts
immediately.

## Empirical calibration (prototype on real `.xst` runs)

The criterion was prototyped against archived pestifer NPT `.xst` files (density taken as
`1/volume`, since the drift/noise metrics are all *fractional* — no PSF mass needed):

- **Soluble target (2bgj, a small solvated protein).** Density jumps 0.906 → ~1.00 within ~1400
  steps of the barostat contracting the box, then just **fluctuates** around ~1.005 for the
  remaining ~17 000 steps of the old ladder — i.e. it is equilibrated by ~step 3 000, and the ladder
  wastes the rest. This is the premise, confirmed. The trend-vs-noise test stops it at **~7–9k
  steps** — versus the ladder's ~20k — a real ~60 % saving, with the mean well-established.
- **Failure modes it exposed.** (a) `fluct_tol` never lets 2bgj converge (intrinsic fluctuations) —
  removed. (b) A pure drift-threshold stops *both* 2bgj and a still-drifting membrane at ~5k, because
  a short window looks locally flat in both — insufficiently discriminating. The trend-vs-noise test
  separates them (membrane pushed to ~22–25k).
- **Scope note.** The membrane npat run drifts ~6 % monotonically over 256k steps and never truly
  plateaus — anisotropic membrane equilibration is **out of scope** for this task (keep the
  npat/npgt protocol). Even trend-vs-noise will eventually declare such a run "converged" over a
  finite window; the `max_steps` ceiling and the residual-drift report are the backstop.

Working defaults from this calibration: `window ≈ 10000`, `n_blocks ≈ 6`, `k ≈ 1.5–2`,
`burn_in ≈ 2000`, `min_steps ≈ 4000`, `chunk_steps ≈ 1000–2000`, `drift_tol ≈ 2e-3` (secondary
ceiling). To be re-confirmed across more systems in validation.

## Parameters (task spec)

| key            | meaning                                                          | default (calibrated) |
|----------------|------------------------------------------------------------------|----------------------|
| `chunk_steps`  | NPT steps per convergence check                                  | 2000                 |
| `max_steps`    | hard ceiling; stop even if not converged                         | 100000               |
| `min_steps`    | never declare convergence before this many steps                 | 4000                 |
| `burn_in`      | leading steps discarded before assessing the trend               | 2000                 |
| `window`       | trailing steps the trend/noise test is computed over             | 10000                |
| `n_blocks`     | blocks the window is averaged into (for the SEM)                 | 6                    |
| `k`            | converge when `trend < k · SEM` (trend insignificant vs noise)   | 2.0                  |
| `drift_tol`    | secondary ceiling: also require `trend/mean < drift_tol`         | 2e-3                 |
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

Run on systems spanning the size range and compare against a long fixed NPT run:

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
  (triple product) + PSF mass; **trend-vs-noise** convergence (`trend < k·SEM`), *not* a fluctuation
  threshold (validated wrong); density-only for v1; a distinct task (not an `md` `until:` mode) for
  discoverability; seeded / deterministic-per-machine; isotropic **soluble-box** scope (membranes
  keep their npat/npgt protocol); replace the ladder across template + `new-system` + examples; keep
  plain `md` for fixed schedules.
- **Open:** final `k`/`drift_tol`/`window` defaults (calibrated on 2bgj; re-confirm across more
  systems); the convergence-report/plot artifact format; whether very large systems want a
  larger `window`/`n_blocks` (better SEM) automatically.
