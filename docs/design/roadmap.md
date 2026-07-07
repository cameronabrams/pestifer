# Pestifer roadmap

A running checklist of planned/possible future work. Not a commitment or a schedule —
just somewhere to park ideas so they aren't lost. Move items into a design doc under
`docs/design/` when they're ready to be worked on; check them off (or delete) when done.

## Solvation

- [ ] **Non-water slabs in `bilayer_embed.tcl`.** The standalone `solvate` task can now
      tile a non-water `kind: box` entry (v3.1.0); extend the same `-spsf/-spdb/-ws/-ks`
      path to the water slabs `bilayer_embed.tcl` lays down, so membrane builds can use a
      non-water solvent too. (Follow-up from `docs/design/solvent-collection.md`, step 6.)
- [ ] **Ship curated built-in solvent boxes.** Build + install a small set of common
      non-water solvent boxes (e.g. MEOH, ETOH, DMSO) into the `solvent` collection so
      users don't each have to run `make-pdb-collection solvent` first.

## Ligands / force field

- [ ] **Automagic CGenFF ligand topologies.** Auto-generate CGenFF topology + parameters
      for an unknown ligand residue encountered in the psfgen task (protonate → mol2 →
      `cgenff` → parse `.str` → remap atom names → inject into psfgen). Design decisions
      are settled (see the `project-cgenff-ligand-workflow` note); not yet built. Adds a
      `cgenff` runtime dependency (free academic binary from SilcsBio).
- [ ] **Automatic update to a new CHARMM force-field release.** Automate ingesting a newly
      published CHARMM release into the bundled `charmmff` resources and re-deriving
      everything downstream of it — segtype classification, atomselect macros, PDB
      collections — instead of the current manual `modify-package charmmff regenerate-*`
      steps.

## Tooling / packaging

- [ ] **Ledger for `modify-package`.** Keep a record of the modifications `modify-package`
      applies to the package (custom residues, PDB-repository entries, force-field changes),
      so contributions are tracked, reviewable, and reversible rather than only visible as
      individual git commits.

## Ring-piercing

- [x] **Lower the `ring_check.cutoff` schema default** (10.0 → 4.0 Å) toward the
      `PSFRing.pierced_by` gate (3.5 Å). Verified result-identical at 3.5/4.0/10.0 on the
      known-piercing fixtures and the ex17 embedded membrane; ~3–7× faster whole-system scan.
      Cutoff-invariance regression test added. (v3.1.0+, `[Unreleased]`.)
- [ ] **Detect ring-piercings at glycan graft time.** Fold a winding-number piercing test
      (reuse `RingChecker` in `psfutil/psfring.py`) into the graft-time declashing
      (`PsfgenTask.declash` / `declash.tcl`), so a threaded glycan is rotated out before the
      structure is written — making the up-front `ring_check` task unnecessary. (See the
      `project-glycan-graft-piercing-check` note.)

## Ideas / unsorted

- [ ] _(add items here)_
