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
