# Design: physics-based loop modeling for missing internal segments

Status: **planned** — direction fixed (fully offline, in-house, physics-based; no
structure predictor and no network). This doc is the plan; nothing built yet.

## Problem

When a source structure has **missing internal residues** (a disordered/unmodeled
loop between two resolved anchors), pestifer builds the residues from sequence and
then closes the gap by steered MD:

1. **psfgen `guesscoord`** builds the missing residues from CHARMM internal
   coordinates — a geometrically valid but essentially arbitrary chain hanging off
   the last resolved residue, its far end *not* bonded to the downstream anchor.
2. **declash** rotates the loop backbone only enough to remove steric overlaps.
3. **ligate/steer** drags the loop's C-terminus across space onto its partner
   N-terminus and applies a peptide-bond patch (`connect`).

The result is topologically correct but **conformationally unrealistic**: steering
optimizes *gap closure*, not *loop plausibility*, so it yields a taut, stretched path
between anchors. `guesscoord` provides no useful starting shape, and declash + steer
cannot recover one.

`ligate` exists **only** for this: `ligate.do()` is entirely gated on
`has_protein_loops` and does nothing else. So this feature replaces `ligate`'s job,
reusing only its final `connect` (peptide-bond patch).

## Goal / non-goals

- **Goal:** produce *plausible, low-energy, clash-free, properly closed* loop
  conformations for missing internal segments, fully offline and reproducibly.
- For genuinely **disordered** loops there is no single "true" conformation; the goal
  is a physically reasonable, relaxable starting structure (an ensemble is the honest
  representation of the ambiguity).
- **Non-goals (for now):** missing *terminal* tails (no downstream anchor → no closure
  problem; sample + minimize, handled separately later); structure prediction / ML;
  any network dependency.

## Approach

Decouple the two things the current design conflates:

- **Conformation source** — where the loop's *shape* comes from.
- **Closure/refinement** — how that shape is stitched onto both resolved anchors with
  correct peptide geometry.

New pipeline, replacing `guesscoord → declash → steer` for internal loops:

1. **Sample** backbone dihedrals (φ/ψ) for each loop residue from a **torsion
   library** (see below) → locally realistic starting geometry instead of `guesscoord`.
2. **Close** the loop onto the downstream anchor with **CCD** (cyclic coordinate
   descent): iteratively adjust the loop's own backbone torsions to bring its
   C-terminal N/CA/C onto the first resolved residue after the gap. Closure is
   distributed across all loop dihedrals — no single end dragged across space (the
   source of the steered-loop artifact).
3. **Filter & score** — reject clashers (reuse the existing declash clash machinery)
   and Ramachandran outliers; rank survivors.
4. **Refine** the best with the existing `minimize` (+ optional short restrained MD),
   then apply the `connect` peptide-bond patch.

## The torsion library

A protein backbone is, to first order, two dihedrals per residue: φ (about N–CA) and
ψ (about CA–C). Their observed joint distribution — the Ramachandran plot — is highly
non-uniform: a few densely populated basins (α, β, PPII, left-handed α) separated by
sterically forbidden regions. A torsion library tabulates that distribution so we can
**sample realistic φ/ψ** for each loop residue.

Refinements captured as **categories** (separate distributions):

- **Residue type:** `general`, `Gly` (permissive), `Pro` (narrow φ), `pre-Pro` (the
  residue preceding a proline has its own preferences).
- **Coil only:** derive from residues **not** in regular secondary structure — a
  *coil library* — so a loop residue isn't sampled as if mid-helix.
- **(P3) Neighbor dependence:** distribution shifts with flanking residue types
  (cf. Ting/Dunbrack neighbor-dependent Ramachandran distributions).

**Representation:** per category, a 2-D histogram of (φ, ψ) on a grid (e.g. 10°×10° →
36×36 cells) storing a normalized probability per cell. A few KB per category; the
whole library is a small versioned data file shipped in the wheel.

**Sampling:** draw a cell weighted by its probability, then a uniform angle within the
cell → a realistic (φ, ψ). Do this per loop residue for a locally plausible backbone,
then CCD closes it.

**Derivation (offline, self-contained — DECIDED):** pestifer derives its **own** coil
library from a curated set of high-resolution PDB structures (for each coil residue,
record type + φ + ψ; histogram per category). A one-time generation script, versioned
in the package like the other regenerated force-field resources — nothing downloaded
at build time.

## CCD closure

Cyclic coordinate descent (Canutescu & Dunbrack, *Protein Sci.* 2003). The loop's
three C-terminal end atoms **M** = {N, CA, C of the last loop residue's junction} must
land on the fixed **target** atoms **F** on the downstream anchor. Sweep each rotatable
backbone bond from N→C; for each bond axis there is a **closed-form optimal rotation θ***
that minimizes S = Σₖ wₖ‖R(θ)Mₖ − Fₖ‖² (rotation of the moving atoms about that one
axis). Apply θ*, advance to the next bond, repeat sweeps until the end-gap RMSD < tol
or max iterations. Constrain/clamp θ to keep φ/ψ in allowed basins, or apply a
Ramachandran filter after closure.

- **P1** uses CCD (simple, deterministic, good MVP).
- **P3** upgrades to analytic **KIC** (kinematic closure): solve the 6 pivot torsions
  that close the loop exactly while non-pivot torsions are library-sampled — more
  diverse, higher-quality sampling (up to 16 solutions per pivot set).

## Determinism (hard requirement)

pestifer values reproducible builds (release-keyed caches; workflow scripts even ban
`Math.random`). The φ/ψ sampler and any stochastic step **must take an explicit seed**
so a config rebuilds bit-identically. The seed is a task spec field with a fixed
default.

## Scoring (two-tier)

- **Rank the ensemble** with a fast internal score (clash count via existing declash
  machinery + Ramachandran favorability, optionally H-bond/SASA) — cheap, no NAMD.
- **Refine only the best K** with NAMD `minimize` (per-candidate NAMD is too slow to
  run on the whole ensemble).

## Task integration

- New loop-building step **owns the gaps** the molecule already tracks
  (`molecule.write_gaps`, `loop_counts`, `has_protein_loops`) and **absorbs `ligate`'s
  `connect`** peptide-bond patch.
- Reuses the existing backbone-rotation Tcl (`write_protein_loop_lines`, `crotations`)
  as the CCD rotate-about-bond primitive; reuses declash for clash scoring and the
  existing `minimize` for refinement.
- **`steer` migration (DECIDED):** the new builder replaces steering for internal
  loops. Keep `steer` as an **opt-in fallback method** for one release (safety net +
  head-to-head comparison on the benchmark), then remove it.
- Likely shape: a `loops:` mod / method selector in the psfgen task (or a dedicated
  post-psfgen task) with `method: ccd` (default) | `steer` (fallback), an ensemble
  size, a seed, and scoring/refinement knobs. Exact schema TBD.

## Validation (first-class deliverable)

Objective benchmark + regression: take structures with **resolved** loops, delete each
loop, rebuild it, and measure **loop backbone RMSD to native** (use the classic
8-/12-residue loop benchmark sets). This gives a quality number, proves P1 beats
steering before we ship, and becomes a regression test. Report per-length RMSD
distributions.

## Phasing (each shippable)

- **P1 — CCD closer replaces steering.** Analytic Ramachandran basins (zero data) for
  the initial φ/ψ, CCD close, minimize, `connect`. Deterministic; smallest change,
  biggest visible win. `steer` demoted to opt-in fallback.
- **P2 — derived coil torsion library + ensemble ranking.** Ship the generated
  `general/Gly/Pro/pre-Pro` coil library; sample M seeded candidates, CCD-close, filter
  clashes, score, keep best-K, minimize.
- **P3 — KIC + neighbor-dependent library + optional restrained-MD refinement.**
  Remove the `steer` fallback once P1/P2 are validated.

## Decisions

- **Fixed:** offline / in-house / no predictor / no network; derive pestifer's own coil
  library from high-res PDB (self-contained, versioned); new builder absorbs `ligate`
  (`connect` kept, `steer` → one-release fallback → removed); CCD before KIC; internal
  loops only; seeded/reproducible sampling; delete-and-rebuild RMSD validation.
- **Open:** exact schema (dedicated task vs psfgen `loops:` mod); grid resolution &
  category set for the library; the high-res PDB subset + culling criteria for
  derivation; ensemble size / best-K defaults; internal score terms; whether P1 ships
  the ensemble or a single CCD closure.

## Scalability analysis (empirical, 2026-07-17)

Motivated by the observation that many bundled examples carry **long** internal loops,
not the 5-residue case P1 was benchmarked on. Missing-internal-loop survey (REMARK 465,
consecutive missing residues flanked by resolved anchors):

| structure | longest internal loop | others | example |
|-----------|----------------------:|--------|---------|
| 7txd | **63** | 18–21 | 11 |
| 4zxb | 37 | 16, 15, 11 | 14 |
| 4zmj / 4tvp | 21 | | 07 / 08 |
| 7xix | 20 | 11 | 15 |
| 5fkw | 18 | | 18 |
| 8fad | 15 | | 09 |
| 5vn3 | 14 | | 12 / 25 |

(HIV-Env glycoproteins: V1/V2, V3, V4 loops.) So the real workload is 14–63-residue
buried loops, far outside P1's tested regime.

**Method.** Engine-level delete-and-rebuild on fully-resolved myoglobin (1mob, chain A,
no gaps): carve internal windows of length N, rebuild with the loop_ccd engine, measure
closure, RMSD-to-native, and steric validity vs the frozen protein environment. Scripts
under scratch (`scaling.py`, `severity.py`, `gen*.py`).

**Findings.**

1. **CCD closure scales indefinitely — never the bottleneck.** Every length 4→28 closes
   to ~0.1 Å in a handful of iterations; iterations *drop* as loops grow (more DOF).

2. **RMSD-to-native saturates at ~10 Å for N≥12** (4→1.5, 8→5.3, 12→8.8, ≥16→~9–11 Å).
   A single basin sample lands in an arbitrary valid closure. (For these disordered
   loops native is ill-defined anyway; the real bar is steric validity, not RMSD.)

3. **The generator, not CCD, is the bottleneck, and its clashes are topological.**
   Deepest heavy-atom overlaps sit at **~0.35 Å** (atoms superimposed = one strand
   passing *through* another), min non-adjacent Cα down to 2.4 Å. These are **not**
   relaxable by minimization/MD (a local optimizer cannot pull one strand back through
   another). Loop-vs-core clashes dominate and appear from **N≈8 up**, not just long
   loops — because per-residue φ/ψ sampling is **environment-blind** and builds the loop
   straight into the fold. Short surface loops (e.g. the BPTI 24–28 benchmark) escape
   only by luck; note that benchmark checked bond+RMSD but **never clashes**.

4. **Closure and sterics are coupled — no sequential heuristic works** (validated across
   four prototypes):

   | approach | closure | sterics |
   |----------|---------|---------|
   | naive sample + full CCD | ✅ exact | ❌ catastrophic, topological |
   | env-aware growth + full CCD | ✅ | ❌ CCD re-swings loop back through core |
   | env-aware growth + clash-aware tail CCD | ❌ 2–11 Å | ✅ soft (~1.3 Å) |
   | land-close growth + 3-point CCD tail | ✅ 0.5 Å (N≥12) | ❌ single arc sweeps core |

   Growth **cannot land the end near the anchor** (median end-gap 7–14 Å): dodging the
   core pushes the end away, and a short pivot cannot bridge that. Exact closure of that
   gap with 3 residues forces a large arc; a single-solution solver (CCD) finds one arc,
   which sweeps the crowded core → many clashes.

**Conclusion — the required architecture.** Exact closure alone is insufficient; the
missing mechanism is **solution diversity at the closure step**. This is precisely
**KIC** (Coutsias–Seok–Dill 2004): analytic closure of a 3-residue pivot yields up to
**16 discrete** closures, giving many distinct arcs to filter for a clash-free one — the
one thing CCD (single local solution) cannot provide. The working generator is therefore
the *combination*:

> **environment-aware self-avoiding growth** (bulk, sterics-only) **+ KIC exact
> multi-solution closure** (3-residue pivot) **+ clash-filtered ensemble** (over
> grown-starts × ≤16 KIC roots) **+ best-K + minimize.**

No single piece suffices; long/buried loops need all of them, at ensemble scale.

**Revised phasing implication.** KIC moves from "P3, nice-to-have" to **required** for
the medium/long (14–63 res) example loops. P1's single-sample CCD is honestly a
**short, solvent-exposed loop** tool (≲8 res); it must not silently emit topologically
broken long loops. Near-term: (a) add a clash/threading diagnostic so the builder and
benchmark fail loud instead of shipping interpenetrated loops; (b) build the native KIC
generator per the architecture above.

## KIC generator: prototyping results (2026-07-17)

Item (a) shipped (`loop_clash_report`/`heavy_env_coords_from_pdb` + `on_clash` knob).
Item (b) was prototyped offline (scratch `kic_*.py`, `distpivot*.py`) with these results:

- **Forward kinematics (NeRF) + internal-coordinate extraction: exact** — reconstructs a
  real backbone to 1.6e-14 Å. Solid substrate.
- **Multi-start closure solver (numpy/scipy `least_squares`, no analytic polynomial):
  correct** — recovers native torsions to 0.000°, all closures exact to ~1e-13 Å, and
  the distinct-solution count is stable across start budgets (2–4 real closures for
  compact tripeptides; not silently missing roots). This is the "robust numerical,
  revisit the Coutsias polynomial later" route.
- **Pivot placement matters decisively.** Reserving the *last 3 consecutive* residues as
  the pivot fails on long loops: the end-gap after env-aware growth is 7–17 Å, more than
  3 consecutive residues can reach → **zero closures**. The canonical KIC choice —
  **3 pivots distributed across the loop** (¼, ½, ¾) — gives each pivot a long lever
  (the whole downstream loop swings), so even a 15 Å gap closes with small torsion
  changes. Distributed pivots restored closure (24–48 solutions/seed) and cut clashes
  **~10×** vs the single-arc CCD (N=20: best ~6 contacts vs 52–109).
- **But not clash-free for crowded loops with a small ensemble.** Delete-and-rebuild on
  myoglobin (all-α, compact; even its most-exposed 20-mer has ~100 env atoms within 5 Å
  of the loop path) still left **deep** clashes in the best of 8 grown-starts × 16
  closures (worst heavy-atom overlap ~0.5 Å, min non-adjacent Cα ~2.5 Å = residual
  threading). Clash-free yield 0 at N=20/28.

**Assessment.** Distributed-pivot KIC is confirmed as the *correct closure mechanism*
(reachability solved, large clash reduction, exact + diverse). Reaching *clash-free* for
crowded medium/long loops needs more than growth + closure + a small ensemble: (i) much
larger ensembles, (ii) a better clash-aware growth than the crude single-pass greedy
walk (proper backtracking self-avoiding growth, biased sampling), and/or (iii) an
energy-based refinement stage. That is a genuine loop-modeler effort with uncertain
payoff on the hardest (buried) loops — which are hard for dedicated tools too. Myoglobin
is a demanding proxy (no truly exposed long loops); the real HIV-Env V-loops are more
solvent-exposed and would likely fare better, but that is untested.

**Open decision.** Whether to (1) invest in the full sampling+refinement modeler,
(2) ship distributed-pivot KIC as a *better-but-not-guaranteed* generator feeding NAMD
minimize (accepting it won't fix deep clashes), or (3) treat long/buried loops as an
external-import problem (Modeller/Rosetta/AF → pestifer stitches + minimizes) and keep
the shipped diagnostic as the guardrail.

## Scope correction (2026-07-17) — this was the wrong problem

The analyses above chase *native reconstruction* of resolved loops (delete-and-rebuild
RMSD-to-native, clash-free re-routing through a packed fold). **That is not the actual
requirement.** In the bundled examples the unresolved internal gaps are, without
exception, **floppy solvent-exposed surface loops** — they are unresolved precisely
*because* they interconvert among many conformations. There is no unique native to
recover, and no need to build a loop-modeling algorithm (growth, KIC, distributed pivots)
at all.

The real requirement is modest: replace the legacy *grow-a-straight-chain-then-steer-the-
end-to-the-anchor* closure with something **more physically defensible** whose one
essential virtue is preserved — **it introduces no new steric clashes**. A straight chain
projects cleanly into solvent; the improvement is a *realistic* backbone that still does.

**Adopted solution (simple, shipped):**

1. Seed each loop residue from **Ramachandran basins** (realistic φ/ψ, not a straight
   line) and close by **CCD** (distributes the small closure adjustment across the loop's
   own torsions — no steered-MD artifact).
2. **Clash-filtered ensemble** — close a handful of independent seeds per gap and keep the
   one that adds the fewest clashes (via `loop_clash_report`). Surface loops have a solvent
   halfspace, so a few seeds reliably include a clean one. This restores the old approach's
   no-new-clashes guarantee with a far better backbone. (`ligate.ccd.ensemble`, default 10.)
3. **Light minimization** of the closed loop (rest of the structure fixed) to relax
   residual torsional strain — a defensible finish, replacing steered MD.

**Retired:** the KIC generator, environment-aware growth, and distributed-pivot closure.
The validated findings above are kept as an honest record of *why* that machinery is
unnecessary here (it targets buried-loop native reconstruction, which the real workload
never asks for). The clash/threading diagnostic stays as a guardrail. If a genuinely
buried, non-floppy long loop ever arises, revisit the retired analysis.
