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

## Open questions for the user

- **Repulsive potential:** hard-sphere with CHARMM vdW radii (simplest, reproducible) vs soft WCA-style
  (better acceptance in dense confinement). Leaning hard-sphere for v1.
- **Cylinder diameter default** and whether to expose it per-species (mixtures with very different
  intrinsic areas may want per-species diameters).
- Whether the headgroup should be confined too, or only the tails.
- Validation target: a bilayer-extracted DMPC conformer footprint/order-parameter reference to check
  against (ask the user; they know the lipid physics and their DB).

## Key code locations

- Packer + conformer selection: `pestifer/molecule/bilayer.py` (`write_grid_pdb`, `bag_of` ~L544,
  `_load_conformer` ~L489, `conf` defaults L159/343/387/549).
- Conformer generation + cache: `pestifer/charmmff/autocache.py` (`ensure_lipid_conformer` ~L185,
  `do_resi`), cache root `~/.pestifer/pdbrepository`, base repo tarball
  `pestifer/resources/charmmff/feb26/pdbrepository/lipid.tgz`.
- Repository entry API: `pestifer/charmmff/pdbrepository.py` (`get_pdb`, `get_conformer_data`,
  `info['conformers']`).
- User's lipid DB (source of the conformers): `/mnt/storage1/cfa/research/lipids/db/`.
- Validation harness: build ex16 (`pestifer build-example 16`) under `~/devtests/pestifer/`; area/APL
  from the `.xst` via `pestifer.util.density_convergence` (`xst_cell_areas`, `membrane_leaflet_geometry`).
