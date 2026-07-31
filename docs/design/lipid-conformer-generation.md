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

## Open questions for the user (blocking final Lever 2a tuning)

- **Cylinder sizing.** The confinement area is *not* the packed APL. What target should map to the
  confinement radius — a measured single-conformer footprint from a real bilayer (per species), or an
  APL×inflation-factor? Per-species (mixtures) or global? The `cylinder_apl` knob and its `sqrt(area/π)`
  mapping are in place; only the calibration is missing.
- **Validation target + reference data.** Preferred fluid-likeness metric — Scd order parameters vs
  mean tail length/footprint — and a reference to hit (e.g. DMPC conformers extracted from a fluid
  bilayer in the user's DB at `/mnt/storage1/cfa/research/lipids/db/`). This closes phasing step 3.
- **Sampling ergodicity.** Accept the current small-move MC with more steps, or add a
  configurational-bias / tail-regrowth move for the dense regime? (Affects `mc_n_equil`/`mc_n_decorr`
  defaults and possibly a new move type in `run_mc`.)
- (Resolved on paper 2026-07-31: hard-sphere potential; tails-only confinement — both implemented.)

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
