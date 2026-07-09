# Pestifer roadmap

A running checklist of planned/possible future work. Not a commitment or a schedule —
just somewhere to park ideas so they aren't lost. Move items into a design doc under
`docs/design/` when they're ready to be worked on; check them off (or delete) when done.

## Solvation

- [ ] **Non-water slabs in `bilayer_embed.tcl`.** The standalone `solvate` task can now
      tile a non-water `kind: box` entry (v3.1.0); extend the same `-spsf/-spdb/-ws/-ks`
      path to the water slabs `bilayer_embed.tcl` lays down, so membrane builds can use a
      non-water solvent too. (Follow-up from `docs/design/solvent-collection.md`, step 6.)
- [x] **Ship curated built-in solvent boxes.** MEOH, ETOH, and DMSO boxes (nmol=216,
      production NPT) are installed in the built-in `solvent` collection (`feb26`), so
      `solvate: {solvent: MEOH}` works without building a box first. (v3.2.0.)

## Membranes

- [x] **Vector-based orientation for membrane embedding.** The `manipulate` task's `transrot`
      mod gained an `ALIGN` movetype that rotates a fragment so a source vector is carried onto
      a target vector (each a literal 3-vector or an atomselection pair), via the minimal
      roll-free rotation about the fragment COM. So a protein can be oriented before embedding by
      aligning a chosen spanning axis onto the membrane normal `[0,0,1]`. (v3.3.0.)
  - [ ] **Follow-up: wire it into the embed step.** Today the ALIGN orientation is a manual
        pre-step in a `manipulate` task; expose it directly as an orientation option in
        `make_membrane_system`'s embed step so a single task can both orient and embed.

## Resources / on-demand generation

- [ ] **Generate-on-miss into a user cache (`~/.pestifer/`).** When a build references a
      residue that is defined in a CHARMM topology/stream file but has **no PDB-repository
      entry**, generate the needed coordinates on the fly and cache them per-user instead of
      hard-erroring (today: `PestiferBuildError('Cannot find {l} in PDB repository')` in
      `bilayer.py`, and the `_solvent_box_args` error in the solvate task). Infrastructure
      mostly exists: `~/.pestifer/` is already the user-resource home (`~/.pestifer/toppar`),
      `charmmffcontent` already accepts `user_pdbrepository_paths` (from
      `charmmff.pdbcollections`), and the generators exist (`do_resi`, `make_solvent_box` +
      the obabel `--gen3d` fallback).
  - The artifact *kind* is driven by the **consumer**, not the species:
    grid membrane packer → single-molecule conformers (`kind: molecule`, via `do_resi`);
    `solvate` task → pre-equilibrated box (`kind: box`, via `make_solvent_box`). Wire the
    solvent-box miss first (that path is freshly built + hardened), then the lipid-conformer miss.
  - Decisions to settle: cache at `~/.pestifer/pdbrepository/<release>/` (release-keyed —
    params are release-specific), auto-registered as a user collection; treat it like a
    **compile-cache** (first use slow with *loud* logging of the one-time cost, cached after);
    an **auto-gen quality tier** (shorter equilibration than the curated shipped boxes); a
    **config toggle** to opt out (reproducibility-strict / offline users prefer the explicit
    error); and an **atomic write / lock** so concurrent builds don't race the tarball.
  - Complements — does not replace — the contribute flow (`make-pdb-collection` +
    `modify-package pdb-repo add-entry`), which is for getting entries into the *shipped*
    package; this is the per-user convenience path.

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

- [x] **Ledger for `modify-package`.** Append-only record at
      `pestifer/resources/modifications.jsonl`; every mutating command records an entry
      (committed alongside the change). `ledger show` lists them; `ledger revert <id>`
      git-reverts a recorded modification on a fresh branch and curates the ledger.
      (v3.2.0.)

## Ring-piercing

- [x] **Lower the `ring_check.cutoff` schema default** (10.0 → 4.0 Å) toward the
      `PSFRing.pierced_by` gate (3.5 Å). Verified result-identical at 3.5/4.0/10.0 on the
      known-piercing fixtures and the ex17 embedded membrane; ~3–7× faster whole-system scan.
      Cutoff-invariance regression test added. (v3.2.0.)
- [ ] **Detect ring-piercings at glycan graft time.** Fold a winding-number piercing test
      (reuse `RingChecker` in `psfutil/psfring.py`) into the graft-time declashing
      (`PsfgenTask.declash` / `declash.tcl`), so a threaded glycan is rotated out before the
      structure is written — making the up-front `ring_check` task unnecessary. (See the
      `project-glycan-graft-piercing-check` note.)

## Ideas / unsorted

- [ ] _(add items here)_
