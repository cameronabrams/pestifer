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

**Why the ladder is *short runs*, not one long run — this is the load-bearing reason.** NPT on a
freshly solvated box shrinks the cell (the water was placed at a lower-than-equilibrium density).
NAMD fixes its **spatial decomposition (patch grid) and PME grid at the *start* of each `run`**; as
the cell shrinks *within* a run the patches shrink with it, and once a patch dimension falls below
the cutoff NAMD aborts ("periodic cell has become too small for the patch grid"). Each new `run`
(a NAMD restart) **recomputes** the decomposition and PME grid for the current, smaller cell — so the
ladder's staged short runs restart *before* the cell shrinks past the patch **`margin`** (pestifer
runs with `margin: 4` Å and auto-sized PME via `pmegridspacing: 1.0`). The runs are *progressively
longer* because the box shrinks fast early (frequent restarts needed) and barely at all once near
equilibrium. So the chunking is **mandatory for numerical stability**, not merely a convergence
convenience — and a single long NPT run does not just waste compute, it **crashes**.

On top of that, the ladder is **a heuristic with magic numbers baked into every config** that does
not adapt to system size, and it is copy-pasted into the template, the `new-system` pipeline, and
every bundled example.

## Goal / non-goals

- **Goal:** one `density_equilibrate` task that runs NPT as a series of **stability-bounded short
  restarts** (short enough that the shrinking cell never outruns NAMD's patch grid) and **stops
  itself** when the box density has converged (or a step ceiling is hit) — replacing the ladder
  across the board (template, `new-system` interactive pipeline, bundled examples). It automates
  *both* jobs the ladder does by hand: staging the restarts, and deciding when enough is enough.
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

## Approach — stability-bounded chunks, convergence-driven stop

Two independent jobs, cleanly separated:

- **Chunk length is a *stability* constraint** (not a free knob). Each chunk is one NAMD `run`
  (restart), which recomputes the patch/PME grid. A chunk must be short enough that the cell does not
  shrink past the patch `margin` *within* it — otherwise NAMD crashes. Since the shrink rate is
  largest early and → 0 near equilibrium, chunk length is **adaptive**: short early, longer late
  (exactly the ladder's shape, but derived rather than guessed).
- **Stopping is a *convergence* decision**, checked at each chunk boundary (the criterion below).

```
firsttimestep = 0; series = []
chunk = chunk_min                                 # conservative first chunk (box shrinks fastest now)
while total_steps < max_steps:
    try:
        run NPT for `chunk` steps (continuation; firsttimestep = total_steps)
    except NAMD "cell too small for patch grid":  # REACTIVE net (see below)
        chunk = max(min_retry_chunk, chunk // 2)   # roll back to last state, retry shorter
        continue
    series += density_samples_from_xst(this chunk)
    total_steps += chunk
    # STABILITY: size the next chunk from the shrink rate just observed; clamp growth to
    # chunk_growth x this chunk (AIMD-style re-growth after a halving) and cap at chunk_max
    rate_sized = clamp(safe_fraction * margin / max_recent_shrink_rate_per_step, min_retry_chunk, chunk_max)
    chunk = min(rate_sized, chunk_growth * chunk, chunk_max)
    # STOP: density convergence (see criterion)
    if total_steps >= min_steps and converged(series):
        break
register the final state; write a convergence report + density-vs-time plot
```

Sizing the next chunk from the *observed* shrink rate keeps every run safely inside the `margin` with
no guessing, and naturally lengthens the chunks as the box settles. `margin` itself can be raised for
extra slack (at a modest performance cost) if needed. This reuses the existing substrate directly:

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

### Reactive stability net — predictive sizing is backed by crash-catch-and-rollback

The predictive chunk sizing above is a *conservative estimate*, and it deliberately does not try to
model NAMD's exact abort threshold — which is genuinely hard to pin down. NAMD picks the patch grid at
`run` startup from the cell, the cutoff/pairlist distance, and `margin`; the true shrink a process can
absorb depends on all of these, and on **CPU the effective `margin` can differ from ours** (the
GPU-resident block sets `margin: 4`, but a CPU run may see a much smaller effective slack). So instead
of betting correctness on a formula, we make the predictor's job *safety-optional* and add a reactive
net that turns the one failure we care about into a self-correcting loop:

- **Detect.** After a `namdrun` raises, scan the NAMD log for the specific abort *"Periodic cell has
  become too small for … patch grid"* (`is_patch_grid_crash`, matched on the stable `too small` +
  `patch grid` fragments). Any *other* failure is re-raised unchanged — a shorter chunk does not fix a
  RATTLE failure.
- **Roll back — for free.** A crashed `namdrun` **registers no new `state` artifact**, so the pipeline's
  current state still points at the *previous completed chunk* — a clean, pre-threshold checkpoint. We
  do not touch NAMD's own `.restart` files or kill/rollback a live process: the last good chunk *is* the
  rollback point, and the retry simply continues from it.
- **Retry shorter.** Halve the chunk (down to a hard floor `min_retry_chunk`, independent of
  `chunk_min`) and rerun. Fewer steps ⇒ less shrink ⇒ we stop before the cell outruns the grid. A retry
  budget (`max_shrink_retries`) bounds this; if even `min_retry_chunk` steps still crash, the box is at
  the grid limit and no shorter chunk helps — we surface a clear error telling the user to raise
  `margin`.
- **Carry the learned size forward (AIMD).** A halved chunk is not snapped back to `chunk_min`; the
  floor for subsequent sizing is `min_retry_chunk`, and growth is capped at `chunk_growth ×` the chunk
  just run. So the scheme *additively-increase / multiplicative-decrease*-es toward the largest safe
  chunk and tracks it upward only as the box actually settles — crucial on a CPU run where the true
  budget is small, so we don't re-crash every single restart.

This is exactly the "monitor-and-preempt" idea, but using **NAMD's own abort as the just-in-time
signal** and **the previous chunk's end-state as the rollback point** — no live `.xst` supervisor
thread, no signal handling, no kill-race on a half-written restart file. The predictive sizing keeps
crashes rare in the common (GPU, `margin: 4`) case; the net makes correctness robust when the estimate
is wrong. (A separate `parse_patch_grid` reads the grid NAMD actually built and logs the approximate
shrink headroom, purely for observability/tuning — nothing depends on that number.)

Why not a live high-frequency `.xst` monitor that stops NAMD *before* the crash? Because a Tcl `run`
loop inside one process does **not** re-decompose the patch grid (that is why `margin` exists), so any
preemption still bottoms out at a process restart — the same boundary chunking already uses — and NAMD
offers no clean checkpoint-and-exit-on-signal, so live preemption means kill + rollback to a
`restartFreq` checkpoint (coarser than `xstFreq`, with a corrupt-file race). The reactive net gets the
same safety with far less machinery, and the barostat shrink rate is *monotonically decreasing* (a
relaxation), so the predictor's only blind spot — acceleration mid-chunk — essentially does not occur.

## The convergence criterion (the crux)

**The observable is equilibrated when its *mean* has stopped drifting — not when its fluctuations
are small.** Density fluctuations are *intrinsic to the NPT ensemble* and do **not** shrink as the
system equilibrates (a small box legitimately fluctuates ±0.3–0.5 %). So the criterion tests the
**trend against the statistical noise**, i.e. whether the systematic drift across the window is
smaller than the uncertainty in the mean.

*(This corrects an earlier version of this doc that also required the fluctuation amplitude itself
to fall below a threshold — validated wrong against real `.xst` runs: an equilibrated soluble box
never satisfies it, because its fluctuations are intrinsic. See "Empirical calibration" below.)*

The design must contend with two facts, both confirmed against real runs (see "Empirical
calibration"):

- The density **fluctuation amplitude is size-dependent**: `σ_ρ/mean ∝ 1/√N_atoms`, measured from
  ~2.4e-3 (a ~10k-water box) down to ~6.5e-4 (a ~180k-water box).
- For a *stationary* (equilibrated) signal, the windowed slope estimate has `trend/SEM ≈ √12 ≈ 3.5`
  **independent of window length** — so a `trend < k·SEM` significance test has a hard floor at
  `k ≈ 3.5` and can never fire for `k ≈ 2`. (This is the mathematical reason every run read
  trend/SEM ~5–6; the significance form is **rejected**.)

Together these mean a **fixed** `(window, drift_tol)` cannot be size-portable: the noise floor of any
drift metric over a fixed window scales as `1/√N_atoms`, so a tolerance tuned on one box is too tight
for a smaller one and too loose for a larger one. The fix is to **make the tolerance a
size-independent practical choice and let the *duration* adapt to the system's noise.**

**Criterion — practical fractional-drift tolerance with an adaptive, precision-targeted window.**

1. **Burn-in + trailing window.** Discard the first `burn_in` steps, then assess only the **trailing
   `window_frac`** (default 0.5) of what remains. *This trailing part is load-bearing and was
   confirmed necessary against a real run* (validation below): a **full-history** window keeps the
   early densification ramp in view forever, so its slope and its block-mean spread never fall — a
   BPTI box that had plateaued by ~30 k steps still read `drift ≈ 5e-3`, `SEM/mean ≈ 8e-4` at 80 k and
   **never converged**. A trailing window lets the ramp age out (once the run exceeds
   `~1/window_frac ×` the transient, the window is all plateau) while still *growing in absolute
   length* as the run lengthens — so the precision gate below keeps tightening.
2. **Precision gate (this is what makes it size-independent).** Grow the trailing window (by running
   more chunks) until the standard error of the windowed mean density is well below the tolerance:
   `SEM/mean < drift_tol / p` (e.g. `p = 3`). `SEM = std(d̄_k)/√n_blocks` from `n_blocks`
   block means `d̄_k`. A small/noisy box simply accumulates a **longer** window automatically; a big
   quiet box reaches the precision gate quickly. Required independent samples scale as
   `n ≳ 12·(σ_ρ/mean)² / drift_tol²`, i.e. `∝ 1/N_atoms` — the duration carries the size dependence,
   not the tolerance.
3. **Drift.** Once the mean is precise enough to *resolve* it, compute the fractional drift
   `drift = |slope|·(t_last − t_first)/mean` from the line fit to `(t_k, d̄_k)`.
4. **Converged** when `drift < drift_tol` for **`n_consecutive` successive checks** (hysteresis),
   with `total_steps ≥ min_steps`.

`drift_tol` is now a **pure adequacy choice** — "the box density has stopped changing by more than X%
over the averaging window," X chosen for the purpose (a starting point for production MD; ~0.1 % is a
sensible target), *not* a measured noise floor. The precision gate guarantees the tolerance is
resolvable above noise for any system size; the hysteresis defends against the noisy, non-monotonic
drift metric (2bgj dips to 6.6e-4 then back to 4.5e-3 near the plateau, so a single passing check
false-triggers — the k~2 failure).

If `max_steps` is reached first, the task **stops anyway** with a loud warning reporting the residual
`drift` and whether the precision gate was even met (a large residual, or an unmet gate, flags a
system that never settled). A NaN / box-blowup aborts immediately.

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
- **Density fluctuation is size-dependent (measured).** Detrended `σ_ρ/mean` over an 8k-step tail:
  2bgj (~10.5k waters) 2.4e-3, ex20 (~32k) 1.3e-3, 6cm3 (~126k) 7.8e-4, GroEL/1aon (~180k) 6.5e-4 —
  a clean `1/√N_atoms` (2bgj/GroEL ratio 3.7 vs √17 = 4.1). This is *why* no fixed-window tolerance is
  size-portable, and why the criterion uses a **precision gate** so the duration (not the tolerance)
  carries the size dependence.
- **A single long (staged) run does not "pin" a universal tolerance.** At a true plateau the mean is
  constant, so any measured windowed drift is *noise* whose floor scales as `1/√N_atoms` (and with
  window) — one system's plateau residual is not a transferable number. (And there is no *non*-staged
  long run to appeal to anyway — a monolithic NPT `run` crashes on the shrinking box; the reference is
  always a chunked run, see Validation.) The tolerance is therefore an **engineering adequacy choice**
  (~0.1 % density stability), and a reference run's job is to *validate the behavior* — that the
  adaptive scheme fires at the equilibrium density and terminates in reasonable time — not to measure
  the tolerance. The `1/√N` characterization above already tells
  us how the required duration scales, so one or two runs at different sizes suffice to confirm the
  precision gate does its job; more are not needed.
- **Scope note.** The membrane npat run drifts ~6 % monotonically over 256k steps and never
  plateaus — anisotropic membrane equilibration is **out of scope** (keep the npat/npgt protocol);
  the `max_steps` ceiling + residual-drift report are the backstop for a system that never settles.

Defaults (form settled; `drift_tol` an adequacy choice, not a measured floor): `drift_tol ≈ 1e-3`
(≈0.1 % density stability), precision-gate `p ≈ 3`, `n_blocks ≈ 6`, `burn_in ≈ 2000`,
`min_steps ≈ 4000`, `n_consecutive ≈ 3`, `window_frac ≈ 0.5` (trailing half). To be sanity-checked
(fires at equilibrium density, sane runtime) on a small and a large system in validation — not
re-calibrated.

### Validation run — BPTI, done (implementation P1)

The task was run end-to-end on real NAMD (BPTI in a ~15 k-atom TIP3 box, CPU). Findings, all folded
back into the code above:

- **The machinery works.** Adaptive chunking `500→700→1050→1570→2350→3520→5000` (grown ≤`chunk_growth`×
  per step, capped at `chunk_max`), no crashes, density from the `.xst` physically sensible
  (`0.98 → 1.03 g/cc`, plateau by ~30 k steps), report + plot artifacts written, continuation
  (`firsttimestep`) threaded across all chunks.
- **Caught a real bug unit tests could not:** NAMD aborts on a `run` whose length is not a whole
  number of `stepsPerCycle`; the AIMD sizing produced 675. The reactive net *correctly declined* it
  (not a patch-grid crash) and surfaced it. Fixed by `quantize_steps` (all chunk lengths → whole
  cycles). *This is why you validate on the real engine.*
- **Caught the full-history-window flaw:** with a full-history window the plateaued box **never
  converged** (ran to the ceiling); the trailing-`window_frac` window (added above) converges with the
  **default** `drift_tol`/`precision_p` at ~step 62 k (density = the plateau value, no premature stop).
  The stop is conservative for so small/noisy a box (its `SEM/mean` bottoms out near the `drift_tol/p`
  gate) — exactly the intended "small boxes run longer"; a large, quiet box will pass the gate far
  sooner. Defaults were *confirmed*, not re-tuned.

## Parameters (task spec)

| key             | meaning                                                          | default             |
|-----------------|------------------------------------------------------------------|---------------------|
| `chunk_min`     | nominal first-chunk length; crash recovery may go below it       | 500                 |
| `chunk_max`     | longest chunk (once the box has settled)                         | 5000                |
| `shrink_safety` | size each chunk so projected per-dim cell shrink < this·`margin` | 0.5                 |
| `margin`        | NAMD patch margin (Å) — the shrink slack per run; from md specs  | 4                   |
| `chunk_growth`  | max factor a chunk may grow over the previous (AIMD re-growth)   | 1.5                 |
| `max_shrink_retries` | patch-grid crash retries per boundary before giving up      | 5                   |
| `max_steps`     | hard ceiling; stop even if not converged                         | 100000              |
| `min_steps`     | never declare convergence before this many steps                 | 4000                |
| `burn_in`       | leading steps discarded before assessing the trend               | 2000                |
| `window_frac`   | assess the trailing this-fraction of post-burn-in samples (ages out the ramp) | 0.5    |
| `n_blocks`      | blocks the window is averaged into                               | 6                   |
| `precision_p`   | precision gate: require `SEM/mean < drift_tol / precision_p`      | 3                   |
| `n_consecutive` | successive passing checks before stopping (hysteresis)           | 3                   |
| `drift_tol`     | converge when fractional drift `< drift_tol` (adequacy choice)   | 1e-3 (~0.1 %)       |
| `seed`          | NAMD RNG seed (reproducible stop on a fixed machine)             | fixed default       |

`chunk_min`/`chunk_max`/`shrink_safety`/`margin` govern the *stability*-bounded chunk length; the
rest govern the *convergence* stop. Plus the usual NPT knobs forwarded to `namdrun` (`temperature`,
`pressure`, `timestep`, output frequencies).

## Task integration

- New task, YAML header `density_equilibrate`, subclassing `MDTask` (it *is* an MD run, with a
  convergence loop around it). `pipeline_contract`: `requires=(STATE,)`, `provides=(STATE,)`.
- `do()` implements the chunk loop above, calling the existing `namdrun(ensemble='NPT', ...)` per
  chunk with the running `firsttimestep`, accumulating density from each chunk's `.xst`, and wrapping
  each `namdrun` in the reactive patch-grid crash-catch/rollback described above.
- Registers the final `state` fileset (as any MD task does) and two artifacts: a **convergence
  report** (per-chunk drift/fluct, the stop reason) and a **density-vs-time plot** (reuse the
  `mdplot` density machinery).
- All NAMD-free logic (density from `.xst`, chunk sizing, crash detection, patch-grid parsing, the
  convergence monitor) lives in `util/density_convergence.py` and is unit-tested in isolation.
- Schema: a `density_equilibrate` task entry mirroring `md`'s NPT options plus the table above.

## Reproducibility

Pestifer values reproducible builds. With an explicit `seed`, NPT is deterministic on a **fixed
machine/NAMD build**, so the adaptive stop point is itself reproducible — the stop step is part of
the (seeded) trajectory. It is **not** bit-identical across platforms (GPU vs CPU, NAMD versions):
the trajectory, and hence the exact stop step, may differ. The *result* is an equilibrated box
either way, but — unlike the offline numpy work (CCD, declash) — this task is not cross-platform
bit-reproducible, by nature of wrapping NAMD. Documented as such.

## Validation (first-class deliverable)

**The reference plateau must itself come from a *chunked* run — there is no "long fixed NPT run" to
compare against.** The very premise of this task (see Problem) is that a single monolithic `run` on a
freshly solvated, shrinking box **crashes** once the cell outruns the patch grid. So the reference
density is *not* a long fixed run; it is a **long *staged* run** — either the old hand-written ladder
carried out to a large step count, or `density_equilibrate` itself with the stop disabled (set
`drift_tol: 0` or an unreachable tolerance so it runs to `max_steps`). Both restart between chunks and
are numerically stable; only the *stopping* differs. The plateau is read off the tail of that staged
reference.

**Validate the *behavior*, don't "calibrate the tolerance."** `drift_tol` is an adequacy choice
(~0.1 %), not a number to be measured (the plateau noise floor is size-dependent — a single run
can't pin it; the precision gate makes the tolerance resolvable regardless). What the validation
*must* confirm — on a **small** and a **large** solvated box (the two ends of the `1/√N` fluctuation
range) — is that the adaptive scheme (i) stops at a density equal to the staged reference's plateau
(not premature), (ii) reaches the precision gate and terminates in a sane number of steps for both
sizes, and (iii) small boxes correctly run longer than large ones. Systems to run:

- a small soluble protein (fast box, should stop early),
- a large glycoprotein trimer (e.g. an HIV-Env example),
- a membrane patch (anisotropic box; confirm the orthorhombic-volume density is sensible, or use
  the full triple-product volume).

Concretely: for each system, run `density_equilibrate` with the stop **disabled** (unreachable
`drift_tol`) to `max_steps` to establish the staged-reference plateau; then run it again with the
production defaults and confirm (a) the auto-stop density matches the reference plateau within
tolerance, (b) the stop step is sensible (not premature, not the ceiling), (c) determinism on a fixed
machine (same seed → same stop). The two runs share the identical chunked trajectory up to the stop
step, so this is a clean A/B on the *stopping rule alone*. Tune the default tolerances/window from
these runs.

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

- **Fixed:** the chunking is **mandatory for numerical stability** (each chunk is a NAMD restart that
  recomputes the patch/PME grid before the shrinking cell outruns `margin`), with chunk length
  **stability-bounded and adaptive** (short early, longer late) from the observed shrink rate — *not*
  merely a convergence convenience; **predictive sizing is a conservative estimate backed by a
  reactive crash-catch-and-rollback** (detect NAMD's patch-grid abort, roll back to the previous
  chunk's registered state — free, since a crashed run registers none — and retry shorter with
  AIMD re-growth), so correctness does *not* depend on modeling NAMD's exact threshold and a live
  `.xst`-monitor/kill supervisor is unnecessary; convergence-or-ceiling stopping; density from the
  `.xst` cell volume (triple product) + PSF mass; convergence by a **practical fractional-drift tolerance made
  size-independent by a precision gate that adapts the window/duration**, with `n_consecutive`
  hysteresis — *not* a fluctuation threshold and *not* a `trend/SEM` significance test (both
  prototyped and rejected: the fluctuation is intrinsic; trend/SEM has a `√12≈3.5` floor at plateau).
  `drift_tol` is an adequacy choice (~0.1 %), not a measured noise floor. Density-only for v1; a
  distinct task; seeded / deterministic-per-machine; isotropic **soluble-box** scope (membranes keep
  npat/npgt); replace the ladder across template + `new-system` + examples; keep plain `md` for fixed
  schedules.
- **Open:** the exact adequacy value for `drift_tol` (~0.1 % is the proposal) and `precision_p`;
  validation that the precision gate fires at the equilibrium density on a small and a large box; the
  convergence-report/plot artifact format.
