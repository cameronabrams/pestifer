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

Density fluctuates; a naive "last frame ≈ previous frame" test false-triggers on a lucky flat
window and, conversely, never converges if it demands raw-frame stability. The criterion must test
the **trend**, not instantaneous values.

**Recommended: block-averaged drift + fluctuation test.**

1. **Burn-in.** Discard the first `burn_in` steps (the initial barostat transient), e.g. a fixed
   step count or a fraction of `total_steps`.
2. **Window.** Consider the trailing `window` steps (or all post-burn-in samples if shorter).
3. **Block-average.** Split the window into `n_blocks` contiguous blocks; take each block's mean
   density `d̄_k` at its mid-time `t_k`. Block averaging decorrelates fluctuations so the trend is
   assessed on ~independent points.
4. **Drift.** Least-squares-fit a line to `(t_k, d̄_k)`; let `slope` be its slope. Normalize:
   `drift = |slope| · window / mean_density` (the fractional change of density across the window
   implied by the trend).
5. **Fluctuation.** Let `fluct = std(d̄_k) / mean_density` (relative spread of the block means).
6. **Converged** when both `drift < drift_tol` **and** `fluct < fluct_tol`, provided
   `total_steps ≥ min_steps`.

Testing *both* a small trend and a small spread avoids the two failure modes: a slow monotonic
creep (small spread, non-zero slope) is caught by `drift`; a noisy plateau (zero slope, large
spread) is caught by `fluct`.

*(A simpler fallback criterion — the relative change of the running mean between two consecutive
windows below `rel_tol` — is easy to reason about and can be offered as an option, but the
drift+fluctuation test is more robust to a slow monotonic creep and is the recommended default.)*

If `max_steps` is reached without convergence, the task **stops anyway** with a loud warning
(reporting the final `drift`/`fluct` so the user can judge), rather than running forever. A NaN /
box-blowup in the NAMD output aborts immediately.

## Parameters (task spec)

| key            | meaning                                                        | default (proposed) |
|----------------|---------------------------------------------------------------|--------------------|
| `chunk_steps`  | NPT steps per convergence check                               | 2000               |
| `max_steps`    | hard ceiling; stop even if not converged                      | 100000             |
| `min_steps`    | never declare convergence before this many steps              | 10000              |
| `burn_in`      | leading steps discarded before assessing the trend            | 5000               |
| `window`       | trailing steps the drift/fluctuation test is computed over    | 10000              |
| `n_blocks`     | blocks the window is averaged into                            | 5                  |
| `drift_tol`    | max fractional density trend across the window                | 1e-3               |
| `fluct_tol`    | max relative spread of block means                            | 1e-3               |
| `seed`         | NAMD RNG seed (reproducible stop on a fixed machine)          | fixed default      |

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

- **Fixed:** chunked NPT with convergence-or-ceiling stopping; density from the `.xst` cell volume +
  PSF mass; density-only for v1; a distinct task (not an `md` `until:` mode) for discoverability;
  seeded / deterministic-per-machine; replace the ladder across template + `new-system` + examples;
  keep plain `md` for fixed schedules.
- **Open:** default tolerances and window/chunk sizes (set from validation); the block-drift test
  vs. the simpler running-mean-delta fallback (lean block-drift); the convergence-report/plot
  artifact format; non-orthorhombic volume (triple product) handling for membranes.
