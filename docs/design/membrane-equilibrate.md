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

## Two-stage protocol (2026-07-31) — constant-area density settle, then tensionless area

**The single-stage tensionless run co-compacts the area even on a *correctly built* bilayer.** The
2026-07-30 flow below pre-equilibrates the pristine quilt with `membrane_equilibrate`, expecting the
area to relax to the fluid DMPC value (~60 Å²). Re-running ex16 with the fluid MC conformer start
showed it does the opposite:

- Build was **correct**: box sized to protein + margin (117.9 × 118.6 Å²), 233 lipids/leaflet → naive
  APL = 13984/233 = **60.0**, exactly the SAPL 60 target (verified: 466 DMPC = 466/118 = 466 residues,
  clean 233/233 split about the midplane).
- The single tensionless piston then **monotonically compacted** it: APL 60.0 (build) → 56.4 (chunk 1)
  → **50.7 and still shrinking** at 124 ps, with density overshooting to 1.022 g/cc and still rising.
  P–P bilayer thickness measured 40.2 Å (gel range; fluid DMPC is ~35–36) — an independent structural
  confirmation of gel-like over-condensation at 310 K, where DMPC should be fluid.

**Root cause.** A freshly gridded bilayer is built at the correct *lateral area* (~SAPL) but is
**under-dense in *z***: the leaflet z-reservation (the isolated melted-conformer z-extent, ~27.7 Å,
taller than the packed leaflet because tails interdigitate at the midplane) over-sizes the box height,
so ρ starts ~0.92. The excess volume physically lives in z, but a tensionless piston with
`useConstantRatio` removes it from the lateral dimensions too (x:y shrink together, z flexes) — the
density-fixing job drags the area down to a metastable gel-like APL, and recovery from there is slow
(high K_A) and fragile (the area monitor can false-converge on the compressed flat spot). Confirmed by
the `.xst`: box 114×115×89 → 108×109×81, both lateral and z shrank ~10%.

**Fix — decouple the two relaxations (standard membrane protocol, cf. CHARMM-GUI NVT→NPAT→NPγT).**
`membrane_equilibrate` now runs two stages by default (`two_stage: true`):

1. **Constant-area density settle** (`useConstantArea yes`, `useConstantRatio no`; `useFlexibleCell`
   stays on): x,y frozen, only z relaxes, so the excess volume leaves **from z alone** and the
   build-correct lateral area is untouched. Gates on **density alone** (area is pinned, not a live
   observable).
2. **Tensionless area relaxation** (`useConstantArea no`, `useConstantRatio yes`): the original NPgT.
   Density is already at target, so only the small lateral move remains — no collapse. Gates on
   **both** density and area (the monitors are `reset()` at the handoff so stage-2 convergence is
   judged on stage-2 samples only, and the area is given `area_min_steps` of *actual* relaxation before
   it may certify).

Build SAPL stays at 60 (it was always correct; 56/50 were compaction, not the packer). The
per-leaflet-geometry `PosixPath` bug (`_read_coor_xyz` called `.lower()` on a `Path`) is fixed so
stage 2 gates on a real protein-corrected APL rather than the `lipids_per_leaflet` fallback. Set
`two_stage: false` to recover the legacy single-stage run (fine for a post-embed stage already near
density equilibrium).

## Build strategy (2026-07-30, final) — pre-equilibrate the membrane, then equilibrate the assembled system

We tried a **minimize-only quilt** (skip pre-equilibrating the pristine bilayer; let a single post-embed
`membrane_equilibrate` do everything) and **rejected it on the evidence.** M5 (below) was implemented
and validated first, so the trial was a *fair* test of the strategy, and it failed to converge:

**The minimize-only trial (ex16, `m5_embed_2026-07-30`).** The gridded bilayer, only minimized, embedded
and M5-solvated, gave an assembled system at **ρ = 0.795 g/cc** — vs **0.96** for the earlier
pre-equilibrated membrane (same atom-number-density; the difference is entirely the un-compacted
lattice, *not* the water — M5's water was measured at bulk density). The post-embed run then had to make
up a ~26 %-volume compaction, and:
- **Density never cleared the drift gate.** It plateaued *physically* (~1.023 g/cc, local rate decayed
  from +0.044 to ~+0.0005 per 10k steps) but the trailing-window drift stayed ~0.009 — the residual slow
  tail of a huge compaction, integrated over a window that grows with the run, never fell below 0.002.
  Headed for the 500k ceiling without converging.
- **Area co-compacted metastably low, then slowly recovered.** It converged *early* (drift < 1e-3 by
  ~100k) but at ~11.8k Å² (corrected APL ~52.8, low for DMPC ~60) because it was dragged down with the
  shrinking box; as density settled it began re-expanding, which *un-converged* it (drift back to ~6e-3).

**So the earlier "area is the slow, unfairly-tested observable" worry is resolved — backwards.** With M5
giving bulk-density water from step 0, the **area is the *fast* observable**; **density** is the
bottleneck, and only because the minimize-only start is so under-dense. Pre-equilibrating the membrane
removes that: it starts the post-embed near equilibrium (~0.96, denser still with M5), no long tail.

**Final flow (both pre- and post-embed equilibration; M5 throughout):**

- **Symmetric** (e.g. ex16): grid membrane → **minimize → NVT → `membrane_equilibrate`** (pre-embed
  quilt: density+area of the pristine bilayer) → embed → **M5-solvate** → **`membrane_equilibrate`**
  (post-embed: the assembled system). Both runs start near-equilibrium and converge.
- **Asymmetric** (e.g. ex17): **calibration patches** get `membrane_equilibrate` to measure preferred
  APL (sizes the grid); the quilt is pre-equilibrated as above; then embed → M5-solvate →
  `membrane_equilibrate`.

**M5 — as built (whole-box solvate + carve), and it works regardless of quilt strategy.** Rather than
tile a pre-equilibrated water box, `bilayer_embed.tcl` now solvates the **whole embed box in one call**
(a thick fill, so its far-z faces are the periodic boundary and there are no thin-slab skins), keeps only
the water in the gap z-regions outside the membrane span, and trims the seam in the existing
overlap-removal pass. Measured on ex16: the added gap water went from **60–100 % of bulk (skinned)** to
**100–120 % (no skins)**. This fixes the *added-water* under-density; the *membrane* under-density is
what pre-equilibration handles.

**APL reporting (2026-07-30).** APL is now measured **per leaflet and protein-corrected**:
`(A − A_protein,leaflet) / N_lipid,leaflet`, with the lipid count taken **per molecule after embedding**
(not the pristine grid count) and the protein cross-section measured **at each leaflet's lipid plane**
(heavy-atom convex hull, bounded to the lipid slab so extramembrane parts in bulk water don't count). The
old `area / lipids_per_leaflet` charged the protein footprint to the lipids (≈ +2–4 Å²/lipid for ex16,
leaflet-dependent) and, if fed the pristine count, hid it by dividing by too many lipids. See
`membrane_leaflet_geometry`. `lipids_per_leaflet` in the spec is now only a fallback.

**Area criterion — measured on a clean pristine bilayer, first cut implemented (2026-07-30).** The
pre-equilibration quilt gave the clean measurement the co-compacted trial couldn't: on the pristine
DMPC bilayer (no protein, no embedding, no compaction) the **density converges fast and rigorously**
(drift < 1e-3, precision met), but the **lateral area is a soft, slow mode** —

- equilibrated `sigma/mean ≈ 0.55%` (4x the density's 0.14%), and
- integrated autocorrelation time `tau_int ≈ 47,000 steps (~94 ps)`, **~100x the density's**.

At that `tau` the trailing window holds only `n_eff ≈ 4` independent area samples, so the density-strict
gates (`drift_tol` 0.002, precision gate 6.7e-4, reliability 6·tau) are **unreachable** — certifying the
area to density precision would need ~1.4 M steps, and the strict criterion simply never converges
(the quilt was headed for the 500k ceiling "waiting on area"). This is *not* a bad start — it is the
area's intrinsic lipid-repacking timescale. (The barostat is not the lever: `tau_area` is ~470x the
200 fs piston period, i.e. lipid-limited, not barostat-limited — see the 100/50 measurement.)

*First-cut area defaults* (schema `area_*`, sized to this measurement; refine as more systems/force
fields are seen): `area_drift_tol 0.012`, `area_precision_p 2.0` (gate 6e-3), `area_autocorr_reliability
4.0`, `area_burn_in 40000`, and a new **`area_min_steps 80000`** floor. The floor is the key piece:
loosening the tolerance alone converges *prematurely* at ~23k steps (the drift metric reads a momentary
flat spot at the bottom of the fast area drop, ~1.9% below the settled area); the `~1-2·tau` floor blocks
that, after which the area converges at ~90-110k with a residual **~0.7%** below the (still slowly
asymptoting) settled value — the honest practical accuracy for a `tau ≈ 47k` mode. These numbers are
tuned to one DMPC bilayer; the principled generalization (scale the floor/burn-in to a robustly measured
`tau_area`, or an equivalence test that adapts to it) is the follow-up.

**Barostat — 100/50, chosen by measurement (2026-07-30).** A controlled experiment (same embedded,
near-equilibrated system; two 40k-step NPgT runs from the *same* restart, differing only in the Langevin
piston) overturned the naive expectation:

| piston period/decay (fs) | density τ (steps) | **area τ (steps)** | means |
|---|---|---|---|
| 200/100 | 1136 | **9144** | ρ 1.0272, area 12041 |
| **100/50** | 1266 | **2109** | ρ 1.0274, area 12078 |

- **Density is *not* piston-limited** (τ ~10x the period already; the faster piston left it unchanged).
- **The area *is* piston-limited** — the faster lateral piston cut its autocorrelation **~4.3x** with the
  equilibrium density and area unchanged (no artifact). The slow area mode was a *sluggish lateral
  piston*, not (only) intrinsic lipid rearrangement.

So `membrane` (NPgT) barostat is now **100/50** (the isotropic-NPT `barostat` block stays 200/100 —
density sees no benefit). Since area convergence cost scales as `sigma^2 * tau`, this is a direct ~4x
speedup on the bottleneck, and it means the "τ ≈ 47,000" pristine-bilayer figure was drift-inflated — the
clean equilibrated area τ is ~9k at 200/100, ~2k at 100/50 — so the area floor above (`area_min_steps`,
`area_burn_in`) was relaxed accordingly. Whether an even faster piston (50/25) helps further or hits the
lipid floor is an open follow-up.

## Integration — one engine, three sites (original M1–M5 framing; see Revised build strategy above)

- **Post-embed (M3, easiest):** replace the hand-written NPgT ladders in examples 16/17 with a single
  `membrane_equilibrate:` task (keep the restrained minimize/NVT warm-up; replace the NPgT tail) — same
  migration pattern as the P2 protein work.
- **Bare-quilt + calibration-patch (M4, deeper):** `make_membrane_system.equilibrate_bilayer` runs its
  pre-embed relaxation through a subcontroller whose task list is
  `[continuation, ...relaxation_protocol stages..., mdplot]`. Since `membrane_equilibrate` is a
  registered task, the subcontroller can run it exactly like the `md` stages. Concrete change:
  1. **Replace the barostatted tail** of the relaxation protocol (the `NPT/NPAT/NPgT` stages) with a
     single `{membrane_equilibrate: {...}}` stage; keep the `minimize`/`NVT` warm-up fixed. The
     hard-coded default protocol (`equilibrate_bilayer`, currently `minimize → NVT → NPT×3 → NPAT×5`)
     becomes `minimize → NVT → membrane_equilibrate`. `_autostage_protocol`'s fixed restart-ramp
     splitting is then superfluous for that tail (`membrane_equilibrate` sizes its own chunks
     adaptively); keep it only for any fixed barostatted stage a user still lists.
  2. **Retire `_area_convergence`** (the crude first-half-vs-second-half area-drift heuristic): the
     "is the APL trustworthy?" flag becomes simply *did `membrane_equilibrate` converge* (vs. hit
     `max_steps`). This needs `membrane_equilibrate` to register a small `converged` artifact that
     `equilibrate_bilayer` reads back (it already imports the subcontroller pipeline).
  3. **The calibration patch is the headline.** `build_patch` → `equilibrate_bilayer` on the per-leaflet
     patches, so this change automatically makes the **preferred APL** — which sizes the whole grid —
     resolve to an honest area convergence instead of a fixed guess + coarse heuristic. Caveat: a small
     calibration patch is *area-fluctuation-limited* (like BPTI density), so it may need **more** steps
     than the old fixed 40k NPgT to converge honestly — that is the real cost of a trustworthy APL;
     keep `max_steps` sane.
  4. The differential-stress **pressureProfile diagnostic** pass stays a fixed `md` stage (it is a
     measurement on the equilibrated cell, not equilibration).
  Schema note: the `relaxation_protocols.patch`/`.quilt` specs are free-form `type: list`, so they
  already accept a `{membrane_equilibrate: {...}}` stage — no schema change. The real touch-point is
  `equilibrate_bilayer`'s per-stage setup loop, which currently keys on `stage['md']` (skipping
  anything else); it must also configure a `membrane_equilibrate` stage (inject `addl_paramfiles`,
  `useflexiblecell`/`useconstantratio`, and pressure-profile keys) the same way it does `md` stages.

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
- **M3 — migrate examples 16/17 post-embed ladders** to `membrane_equilibrate` (in progress: example 16
  migrated, validation build running).
- **M4 — wire the engine into `make_membrane_system` pre-embed** (quilt + calibration patch); retire
  the fixed NPgT tail of `_autostage_protocol` and `_area_convergence`. **Design drafted** (see
  Integration above): replace each relaxation protocol's barostatted tail with a `membrane_equilibrate`
  stage, make the "APL trustworthy" flag = did it converge, and configure the stage in
  `equilibrate_bilayer`'s per-stage loop (which currently only handles `md`). Calibration-patch APL is
  the headline win. Implement after M3 validates.
- **M5 — equilibrated-density embed solvation (follow-up; not a `membrane_equilibrate` change).**
  *Motivation, measured on the M3+M4 ex16 build:* even with a fully pre-equilibrated bilayer (M4
  converged @119k), the **post-embed** `membrane_equilibrate` densifies very slowly — the lateral area
  settles fast but the **density keeps rising** (c_z shrinks at fixed area), i.e. an under-dense body of
  bulk water slowly compacts. Root cause: `bilayer_embed.tcl` fills the gaps above/below the protein
  with VMD `solvate -minmax {slab}` (lines ~277, ~302). VMD's source water box is equilibrated, but
  **clipping a generic box to arbitrary thin slab bounds leaves an under-dense ~3 Å skin on every cut
  face** (worse for small `zdist`), and the water-box periodicity does not lattice-match the periodic
  system — so the fill is under-dense and the post-embed run has to squeeze out ~1.5% along z over
  hundreds of k steps. *Fix:* **replicate the just-converged pre-embed bilayer's own water slab in z**
  to fill the added gaps instead of `solvate -minmax` — that water is already at this system's exact
  equilibrium density *and* laterally lattice-matched/periodic with the box, so it stacks with no
  trimming skins. Mechanically: carve the bilayer water/ion slab, tile it in z across the gap, trim only
  at the far z-face, splice at the membrane interface, autoionize the added slab as before. Localized to
  `bilayer_embed.tcl` (the bilayer water is already loaded there). Payoff: the post-embed
  `membrane_equilibrate` then converges fast (the real pre/post synergy) because the system starts near
  density-equilibrium. Does not block M3/M4 — those tasks are validated; this is an embed-quality
  improvement that makes membranes converge quickly.
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
