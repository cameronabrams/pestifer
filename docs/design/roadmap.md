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
- [ ] **On-demand solvent boxes: pull in the solvent's own parameter file.** On-demand
      generation (v3.5.0/v3.6.0) equilibrates a box with the *default* parameter set, so a CGenFF
      solvent whose parameters live outside that set fails at the NPT step — e.g. acetonitrile
      (`ACN`) dies with `UNABLE TO FIND DIHEDRAL PARAMETERS FOR HGA3 CG331 CG1N1 NG1T1`. Discover
      and add the RESI's required parameter file(s) before equilibrating (the same way
      `make_solvent_box` already appends the defining topology file), so params-incomplete
      solvents build too. (Discovered building Example 24; acetone works, acetonitrile does not.)

## Membranes

- [x] **Vector-based orientation for membrane embedding.** The `manipulate` task's `transrot`
      mod gained an `ALIGN` movetype that rotates a fragment so a source vector is carried onto
      a target vector (each a literal 3-vector or an atomselection pair), via the minimal
      roll-free rotation about the fragment COM. So a protein can be oriented before embedding by
      aligning a chosen spanning axis onto the membrane normal `[0,0,1]`. (v3.3.0.)
  - [x] **Follow-up: wire it into the embed step.** `make_membrane_system`'s embed step now
        orients via the `transrot` `ALIGN` path: a new `embed.orient` spec (`source`/`target`
        vectors) drives it directly, and the existing `z_head_group`/`z_tail_group` shorthand maps
        onto it (`source: [z_tail, z_head], target: [0,0,1]`). Verified coordinate-identical to the
        former `Orient::orient`-based `bilayer_orient` script (RMSD < 0.01 Å on the van3 fixture).
        (v3.4.0.)

## Resources / on-demand generation

- [x] Generate-on-miss into a user cache (`~/.pestifer/`). When a build references a residue that
  is defined in a CHARMM topology/stream file but has **no PDB-repository entry**, generate the
  needed coordinates on the fly and cache them per-user instead of hard-erroring. Cache at
  `~/.pestifer/pdbrepository/<release>/` (release-keyed), auto-registered as a user collection;
  compile-cache semantics (first use slow with loud logging, cached after); `config toggle`
  `charmmff.generate_missing_coordinates` (default on) to opt out; atomic publish + `fcntl` lock
  so concurrent builds don't race. Complements — does not replace — the contribute flow
  (`make-pdb-collection` + `modify-package pdb-repo add-entry`).
  - [x] **Solvent-box miss** — `solvate` with a non-water solvent that has no shipped box now
        builds one with `make_solvent_box` and caches it under
        `~/.pestifer/pdbrepository/<release>/solvent/<RESI>/` (`kind: box`, marked `quality: auto`).
        (v3.5.0.) Since the decided quality tier is the **full** shipped equilibration (50k NPT),
        auto-boxes are same-quality-but-not-hand-curated rather than a lower tier.
  - [x] **Lipid-conformer miss** — the grid membrane packer's
        `PestiferBuildError('Cannot find {l} in PDB repository')` in `bilayer.py` now generates the
        residue's single-molecule conformers (`kind: molecule`, via `do_resi`) on miss and caches
        them under `~/.pestifer/pdbrepository/<release>/lipid/<RESI>/`, sharing the solvent path's
        lock/atomic-publish/isolated-builder machinery. (The artifact *kind* is driven by the
        consumer, not the species: packer → conformer, solvate → box.) (v3.5.0.)

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

## Glycans

- [ ] **Build-or-extend glycans (`glycosylate` mod).** Today glycans are added by *grafting* a
      pre-existing glycan (residue range from a source PDB) onto a target residue — a `grafts`
      mod in the psfgen task. Add a complementary **`glycosylate` psfgen mod** that *builds* a
      glycan from a specified sequence/linkage tree (and/or *extends* an existing one) rather than
      requiring a donor structure, so a user can request, e.g., "an N-linked complex biantennary
      glycan on ASN 234" declaratively.
  - Decisions to settle:
    - **Specification.** How the user names the glycan — a linear/branched sequence of CHARMM
      sugar RESIs plus the linkage patches between them (N-linked via `NGLB`/sequon chemistry;
      O-linked on Ser/Thr), possibly a shorthand for common cores (Man3/Man5, complex
      biantennary, etc.).
    - **Geometry.** Build coordinates from CHARMM internal coordinates (grow residue-by-residue
      off the attachment point, like a psfgen `patch` + `guesscoord` chain) vs. graft a
      pre-built core then extend — the former needs no donor PDB but leans on IC completeness.
    - **Attachment point.** Reuse the graft machinery's target-residue handling (align onto the
      Asn/Ser/Thr, replace, insert the rest into the topology).
    - **Reuse the graft-time hygiene.** Fold in the existing declash and the planned graft-time
      ring-piercing check (see `project-glycan-graft-piercing-check`) so a built glycan is
      declashed and un-pierced before the structure is written.
  - Relationship to grafting: `grafts` stays for "I have this exact glycan structure";
    `glycosylate` is for "build me this glycan" — they can share the target-attachment and
    hygiene code paths.

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
- [x] **Detect ring-piercings at glycan graft time.** `PsfgenTask.declash` now runs a
      winding-number piercing scan (reusing `RingChecker`) after building glycans and clears each
      glycan-bond piercing by rotating the offending glycan sub-branch out of the ring, before the
      structure is written — making the up-front `ring_check` unnecessary in the common case.
      Rotation engine extracted to `psfutil/ring_resolve.py`, shared with `ring_check`. Toggle
      `glycans.declash.check_piercings` (default on); best-effort/non-fatal. Validated on the real
      4zmj model (2 glycan piercings → 0). (v3.6.0.)

## Ideas / unsorted

- [ ] _(add items here)_
