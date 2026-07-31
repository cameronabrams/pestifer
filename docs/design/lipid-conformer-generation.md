# Lipid conformer generation for fluid-bilayer-like grid packing

**Status:** planned (not started). Handoff doc — self-contained; a new agent should be able to
pick this up cold. Written 2026-07-31.

## Why this project exists

The grid membrane packer places lipids that are **over-packed and conformationally frozen**, and the
membrane_equilibrate task then spends a large budget slowly undoing it. The root cause is *how single
lipid conformers are generated and selected*, not the equilibration machinery (which is committed and
validated — see [membrane-equilibrate.md](membrane-equilibrate.md) and commit `feda3740`).

### What we observed (DMPC, ex16 = HIV gp41 MPER-TM in DMPC)

- A freshly-gridded, equilibrated **pristine** DMPC bilayer settles to **APL ≈ 50 Å²/lipid** (area
  ~11,600 Å² / 233 lipids-per-leaflet), well below fluid DMPC (~60 Å²; CHARMM36 ≈ 60–61, experiment
  ≈ 60.6 at 303 K). It is nearly gel-phase packing. Equilibrated system density ~1.026 g/cc, also high.
- The post-embed area drift is **persistently positive and slow** — the membrane keeps trying to
  expand. Measured lateral-area integrated autocorrelation time on the equilibrated bilayer was long
  (thousands of steps; a drift-contaminated pristine estimate read ~47k). It is the slowest mode in the
  whole build and was the reason the post-embed `membrane_equilibrate` needed ~800k steps / ceilinged.

### The mechanism (confirmed in code + data)

The grid packer (`pestifer/molecule/bilayer.py`, `write_grid_pdb` → `bag_of`, ~line 544) loads **one**
conformer per species — `_load_conformer(nm, specs.get('conf', 0))` — and stamps that *same*
conformation onto every lipid, varying only a random in-plane spin. `conf` defaults to `0` in four
places (bilayer.py lines 343, 387, 489, 549) and the conformer specstring defaults to `'0'` (line 159).

The available DMPC conformers (10 of them, cached in `pestifer/resources/charmmff/feb26/pdbrepository/
lipid.tgz`, generated from the user's DB at `/mnt/storage1/cfa/research/lipids/db/data/lipid/DMPC/` via
`ensure_lipid_conformer` → `do_resi`, a "stretch then relax in vacuum" sampler) are **all extended thin
rods**:

| conformer | z-extension (Å) | xy convex-hull footprint (Å²) |
|---|---|---|
| DMPC-00 (**the one the packer uses**) | 28.7 | **34** (thinnest of all) |
| DMPC-01…09 | 26.9–27.9 | 38–50 |
| DMPC-init (pre-sample seed) | 22.1 | 103 (globular/collapsed — opposite extreme) |

So: the packer picks the single *thinnest* conformer; and even sampling across all 10 only reaches a mean
footprint ~41 — the whole ensemble lacks fluid diversity. Thin rods pack trivially to low APL; the slow
positive area drift is then the acyl chains **acquiring intramolecular configurational entropy** (t→g
isomerization, lateral fattening) toward the fluid state — a slow intramolecular relaxation that the
grid never seeded.

**Why vacuum single-molecule sampling can't fix itself:** in vacuum, attractive vdW drives chain
*collapse* into a globule (see DMPC-init, footprint 103); avoiding that by stretching/restraining
overshoots into thin rods (00–09). Neither is the moderately-splayed, gauche-rich state a lipid holds in
a fluid bilayer, because that state is stabilized by the **membrane environment** (lateral packing +
hydrophobic effect), which is absent for an isolated molecule in vacuum. And we can't reliably generate
real bilayers to sample from — especially for arbitrary **mixtures**. This is the chicken-and-egg the
project must break.

## Design (agreed on paper 2026-07-31)

Break the chicken-and-egg with **athermal (repulsive-only) Metropolis MC**: a purely excluded-volume
potential (no attractive vdW) removes the collapse driving force, so chains explore torsional entropy
freely and you get a melt-like ensemble — a good model for the hydrophobic core, where tail statistics
are packing-dominated. Athermal = geometry-only (T→∞); you are sampling *configurational entropy subject
to excluded volume*, exactly the missing quantity. Hard-core/WCA overlap checks are cheap (no forces;
cell lists → O(N)).

Two complementary **levers**:

### Lever 1 — packer draws per-lipid across the ensemble (quick, do first)
`bag_of`/placement in `bilayer.py`: instead of caching one conformer per species, load **all** conformers
and draw one per lipid from the packer's RNG (already seeded — `rng` in `write_grid_pdb`). Small change,
modest effect on its own (mean footprint ~41 vs 34), but necessary so any improved ensemble is actually
used. Decide whether to include conformer 0 / the stretched seed in the draw (probably exclude artificial
stretches; sample the relaxed ensemble).

### Lever 2a — cylinder-confined athermal MC conformer generator (the real fix; lead with this)
Replace the vacuum-MD sampler inside `ensure_lipid_conformer`/`do_resi` with an athermal MC that samples
tail dihedrals under **lateral cylindrical confinement**:
- **Moves:** dihedral pivot moves on acyl-tail C–C torsions only; **bonds and angles held rigid**
  (preserves CHARMM geometry; MC is pure torsion + excluded volume). Glycerol/phosphate/headgroup
  dihedrals fixed or lightly sampled — the tails hold the entropy.
- **Confinement:** a cylinder whose axis is the head→tail reference axis (the packer's existing
  reference atoms define it), diameter a tunable knob defaulting from a target APL:
  `d ≈ 2·√(APL_target/π)` ≈ 8.7 Å for APL 60. It need not be exact — NPgT sets the final APL; the
  cylinder just prevents rod formation and globule collapse.
- **Output:** same as today — a set of conformer PDBs + `info.yaml`, cached per species → **mixtures are
  trivial** (each species sampled independently; packer mixes). Athermal MC naturally yields a
  *distribution* of footprints, so lever 1 then samples real diversity.
- **Cache key:** extend the existing one with cylinder diameter, move count, athermal flag, seed.

### Lever 2b — whole-grid athermal MC pre-relax (optional follow-on)
Athermal t→g MC (plus whole-molecule rotation/translation) on the **assembled** gridded system, where the
confinement is the *real* neighbor environment (captures packing correlations a mean-field cylinder
misses). Not cacheable per-species — it's a packing relaxation, so it could **replace or precede the
current minimize** (double duty: relieve grid clashes + seed chain entropy) before NAMD runs. Try
cylinder-only first; add this only if the packed membrane still starts too tight.

## Suggested phasing

1. **Lever 1** in `bilayer.py` — per-lipid conformer draw. Ship + a unit test (placement uses >1 distinct
   conformer). Low risk.
2. **Lever 2a** — the cylinder athermal-MC generator as a new module, wired into `ensure_lipid_conformer`
   as the sampler (keep the cache/publish flow). Regenerate the DMPC ensemble.
3. **Acceptance test** (closes on the original bug): generate → grid-pack → NPgT-equilibrate ex16's
   pristine bilayer → confirm **(a)** equilibrated APL matches literature (DMPC ~60), **(b)** area drift
   is small from the start (slow mode largely gone). Also check the generated ensemble's **tail order
   parameter / mean tail-length distribution** against fluid-bilayer values — a more direct
   "are these fluid-like?" check than APL, and it tells you if the cylinder diameter needs adjusting.
4. **Lever 2b** only if needed after measuring (3).
5. Re-tune the area convergence defaults (they were sized to the *over-packed* membrane's slow drift; a
   fluid start should let them tighten — see the `area_*` defaults in `schema/base.yaml` and
   membrane-equilibrate.md).

## Implementation status (2026-07-31)

**Lever 1 — DONE, committed** (`f8772a04`). `bag_of` in `write_grid_pdb` caches the whole conformer
ensemble per species and draws one per lipid (and per chamber solvent) from the packer RNG; explicit
per-species `conf` still pins. Fixed a latent key mismatch (`bilayer_specs['conformers']` vs schema
`lipid_conformers`) that had been silently ignoring the pin knob; schema default `'0'`→`''` so unpinned
ensemble draw is the true default. Unit-tested (`test_bilayer_write_grid_pdb_samples_conformers`).

**Lever 2a — engine + builder + wiring DONE, committed** (`92d4953d`, `44fc6029`, `bd68dffa`); **TUNING
+ VALIDATION OPEN.**
- `pestifer/charmmff/athermal_mc.py`: force-field-agnostic MC core (`MoleculeMC`, `run_mc`,
  `build_lipid_mc`, `cylinder_radius_for_apl`, `build_exclusions`, `moving_set`). Dihedral-pivot moves
  (rigid tip-subtree rotation → bonds/angles preserved exactly); hard-sphere overlap (sigma-scaled
  `Rmin/2`, exclude ≤1-4); **monotone-inward** cylinder confinement; membrane-normal axis through the
  tail-bundle centroid. 16 unit tests (synthetic alkane + synthetic diacyl).
- Wired as opt-in `sampler='mc'` through `do_psfgen`/`do_resi`/`ensure_lipid_conformer` — reuses the
  psfgen-init + minimize + orient front matter, skips stretch/sample MD, then MC → `{resid}-NN.pdb`.
  Provenance recorded in `info.yaml['generation']`. **Default stays `'md'`** until validated. Legacy md
  path verified unchanged end-to-end.

### Cylinder sizing — engineering solution (2026-07-31, agreed with user)
No reference bilayer exists to calibrate against (chicken-and-egg: the user trusts only bilayers they
build, and none exist yet). So sizing is an **engineering knob**: confinement cross-section =
`cylinder_inflation × cylinder_apl`, where `cylinder_apl` is the honest target APL (default 60) and
`cylinder_inflation` (default **1.9**) accounts for a single conformer's convex-hull footprint exceeding
the packed APL. **Cheap calibration proxy** (needs no bilayer): tune inflation so the ensemble's mean
footprint ≈ target APL. On DMPC (shipped MC params), inflation 1.9 (area 114 Å²) → mean footprint ~56–61,
~7/10 distinct shapes, vs the old all-rod ensemble (34–50, mean ~41). Committed `3c242eff`.

**Integration validated end-to-end** (`~/devtests/pestifer/mc_validate/`): generate MC DMPC → cache
(`info['generation']` stamped `sampler:mc`) → grid packer draws per-lipid → 128-lipid patch with
footprints ~45–57. The full chain works. (The manually-seeded MC cache entry was removed afterward to
keep the user env on the shipped md conformers; reseed with `ensure_lipid_conformer('DMPC', CC,
sampler='mc')`.)

### What the DMPC end-to-end runs showed (findings that set up the tuning)
Validation harness: `~/devtests/pestifer/mc_validate/gen_dmpc_mc.py` (`do_resi(..., sampler='mc')`,
then per-conformer xy convex-hull footprint + z-extension).
- The mechanism works: geometry preserved, rod melts partially, footprint compresses and **distributes**
  (e.g. 43–64 Å² at cylinder area 100), distinct conformers produced.
- **Cylinder tightness for two-tail lipids is the key knob.** `d ≈ 2√(APL/π)` (R≈4.37 for APL 60) treats
  a lipid footprint as one disk of area APL, but a *single* conformer's convex-hull footprint runs well
  above the packed APL (lipids interdigitate/tessellate in a real bilayer). At R=4.37 the two tails only
  fit as straight, tightly-stacked rods → the cylinder *selects for* the rod and diversity collapses
  (≈7/10 samples identical). Loosening to area ~90–120 restores diversity. Provisional default set to
  `cylinder_apl=100`, flagged for retune.
- **Dense-MC mixing** limits diversity once compacted: at a tight footprint most pivots clash, acceptance
  drops, consecutive decorrelation windows repeat. Small moves (`mc_max_angle=π/6`) + `radius_scale=0.8`
  gave the best acceptance in a sweep; may still need many more steps or smarter moves (see below).
- z-extension only drops modestly (~29.5→~24–26 Å); not yet clearly "fluid" — needs an order-parameter
  target to judge.

## Acceptance test — RUNNING (phase 3, 2026-07-31)

The bare-membrane contract fix (`58b75144`) unblocked a config-driven pristine DMPC bilayer, and
grid-packing *fluid* (fat) conformers surfaced **two independent failure modes**, both now fixed:

1. **Topology corruption → NaN.** The psfgen step loads the bare grid PDB into VMD to split by
   type/leaflet; with bond perception on, two near-coincident lipid atoms get spuriously bonded and
   VMD **merges their two residues**, scrambling the per-residue coords → psfgen guesses a
   box-sized-displaced atom → NAMD diverges to NaN. Fixed by loading with **`autobonds off`**
   (`c6612865`) — connectivity is irrelevant there (it comes from the RESI topologies at psfgen
   time) and `residue` grouping then falls back to resSeq/segname. (Verified: a 0.05 Å
   near-coincidence loads as 5 merged residues with autobonds on, 6 correct with it off.)
2. **Physical overlap → segfault.** Even with correct topology, fat conformers packed at SAPL 60
   have severe steric clashes that **segfault the NAMD minimize** (rc 139). Fixed by the packer
   **re-spin guard** (`c56863db`): re-spin each lipid until its atoms clear placed lipids by ≥0.9 Å.

Both are needed (correctness + physical viability) and both are committed. With them, a pristine
DMPC bilayer (100 lipids/leaflet) grid-packs from MC conformers and runs `membrane_equilibrate`
(NPgT) to completion.

**First result (inflation 1.9):** equilibrated **APL ≈ 67 Å²** with **tiny drift (−0.52%** over the
final half) — vs the old vacuum-rod membrane's gel-like **~50** with a large *persistent* positive
drift. The slow mode is largely gone (the project's core goal). But 1.9 *overshoots* fluid DMPC
(~60–61): the tails are a bit too fat, so the membrane settles high. Re-running at **inflation 1.5**
(footprints ~48 vs ~56) to bring the equilibrated APL down toward 60. Harness:
`~/devtests/pestifer/mc_validate/` (`pristine_dmpc_equil.yaml`, `analyze_apl.py`).

## Remaining work (superseded — see above)

The only real validator (no reference bilayer) is: **generate MC ensemble → grid-pack → NPgT-equilibrate
→ does APL land ~60 with small drift from the start?** Status/obstacles:

- **A protein-free pristine-bilayer build is currently unreachable via a config file.** `make_membrane_
  system`'s contract intends a bare membrane to be an origin (`requires=()` when not embedding), but the
  ycleptic schema populates the `embed` block with defaults even when the user omits it, so
  `bool(specs.get('embed'))` is always true and the pipeline validator demands a prior protein state.
  Either (a) run the full **ex16** (has a protein — the design's stated acceptance test, but expensive:
  membrane_equilibrate ran ~800k steps for the old over-packed case), or (b) fix the contract to detect a
  *meaningful* embed (e.g. presence of `z_ref_group`/user-set keys) so bare builds work — this is a
  genuine latent bug worth a small separate fix.
- To run the test: seed the MC cache (`ensure_lipid_conformer('DMPC', CC, sampler='mc')`) so the packer
  uses MC conformers, build, and read APL/area drift from the `.xst` via
  `pestifer.util.density_convergence` (`xst_cell_areas`, `membrane_leaflet_geometry`).
- If APL still low / drift still large → nudge `cylinder_inflation` up (looser cylinder → fatter
  footprints) and/or improve MC mixing before considering Lever 2b.

## Still-open tuning knobs (revisit after the acceptance test)

- **Sampling ergodicity.** Current small-move MC gives only ~7/10 distinct conformers (dense-regime
  mixing: compact states have low pivot acceptance). "See how far we get without advanced MC first"
  (user, 2026-07-31). If insufficient: more steps, or a configurational-bias / tail-regrowth move in
  `run_mc`.
- **Per-species inflation** for mixtures with very different intrinsic areas (single global factor now).
- (Resolved 2026-07-31: hard-sphere potential; tails-only confinement; APL×inflation sizing — all
  implemented. No reference bilayer exists, so the acceptance test is the calibration authority.)

## Key code locations

- **Athermal-MC engine (new):** `pestifer/charmmff/athermal_mc.py` — `run_mc`, `build_lipid_mc`,
  `cylinder_radius_for_apl`, `MoleculeMC`. Tests: `tests/unit/test_charmmff/test_athermal_mc.py`.
- **MC sampler wiring:** `pestifer/charmmff/make_pdb_collection.py` — `_sample_and_write_mc_conformers`
  and the `sampler='mc'` branch in `do_psfgen`/`do_resi`; MC knobs (`cylinder_apl`, `mc_n_equil`,
  `mc_n_decorr`, `mc_seed`, `mc_max_angle`, `mc_radius_scale`).
- Packer + conformer selection: `pestifer/molecule/bilayer.py` (`write_grid_pdb`, `bag_of` ~L544,
  `_load_conformer` ~L489; per-lipid draw is now default, explicit `conf` pins).
- Conformer generation + cache: `pestifer/charmmff/autocache.py` (`ensure_lipid_conformer`, threads the
  `sampler`/MC knobs, stamps `info['generation']`), cache root `~/.pestifer/pdbrepository`, base repo
  tarball `pestifer/resources/charmmff/feb26/pdbrepository/lipid.tgz`.
- Repository entry API: `pestifer/charmmff/pdbrepository.py` (`get_pdb`, `get_conformer_data`,
  `info['conformers']`).
- User's lipid DB (source of the conformers): `/mnt/storage1/cfa/research/lipids/db/`.
- Validation harness: build ex16 (`pestifer build-example 16`) under `~/devtests/pestifer/`; area/APL
  from the `.xst` via `pestifer.util.density_convergence` (`xst_cell_areas`, `membrane_leaflet_geometry`).
