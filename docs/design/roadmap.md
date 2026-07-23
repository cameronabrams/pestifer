# Pestifer roadmap

A running checklist of planned and possible future work — somewhere to park ideas so they
aren't lost. It is not a commitment or a schedule: items are checked off as they ship, and
what appears here is refined and reprioritized as the project evolves.

## Solvation

- [ ] **Non-water slabs in `bilayer_embed.tcl`.** The standalone `solvate` task can now
      tile a non-water `kind: box` entry (v3.1.0); extend the same `-spsf/-spdb/-ws/-ks`
      path to the water slabs `bilayer_embed.tcl` lays down, so membrane builds can use a
      non-water solvent too. (Follow-up from `docs/design/solvent-collection.md`, step 6.)
- [x] **Ship curated built-in solvent boxes.** MEOH, ETOH, and DMSO boxes (nmol=216,
      production NPT) are installed in the built-in `solvent` collection (`feb26`), so
      `solvate: {solvent: MEOH}` works without building a box first. (v3.2.0.)
- [x] **On-demand solvent boxes: robust parameter provisioning.** `make_solvent_box` now (A) loads a
      residue's companion parameter file when it is defined in a bare `.rtf`, and (B) turns a
      missing-parameter NAMD abort into an actionable error naming the exact residue and atom-type term
      instead of a cryptic mid-equilibration fatal. Genuinely-degenerate torsions (e.g. acetonitrile's
      `HGA3 CG331 CG1N1 NG1T1`, degenerate because the C–C≡N skeleton is linear) are covered by a
      curated, reviewed `k=0` entry in `custom/toppar_pestifer_dihedral_fills.prm` — so `ACN` now boxes
      — while unknown gaps are still surfaced, never silently zero-filled. Design doc:
      `docs/design/on-demand-solvent-parameters.md`. (Discovered building Example 24.)
- [ ] **Solvent blends / mixtures.** `solvate` takes a single `solvent:` species today; add support
      for a **multi-component solvent** at a specified composition — e.g. a water/DMSO cosolvent mix
      or a ternary blend — for mixed-solvent / cosolvent-MD builds. Decisions to settle:
    - **Specification.** How the user names a blend and its ratio — a list of `{solvent, fraction}`
      entries (mole vs. volume vs. number fraction — pick one and be explicit), with a rule for the
      remainder/primary species.
    - **Box generation.** Pack a mixed box directly (extend `make_solvent_box` to place multiple
      species at the target composition), vs. tile/overlay single-species boxes — the former gives a
      properly mixed liquid; reuse the on-demand generation cache, keyed by the (species, ratio) blend.
    - **Parameters + neutralization.** Union the defining topology/parameter files of all components
      (composes with the on-demand-parameter item above), and keep ion neutralization working against
      the mixed box.

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

- [x] **GPU mode detection shouldn't hinge on path inequality.** The `namd.processor-type`
      schema entry (formerly deprecated-and-ignored) is now a functional `auto|cpu|gpu` control
      (default `auto`). `auto` keeps the historical path-inequality detection (+ the
      `_warn_gpu_not_elected` warning); `gpu` **forces GPU mode even in the single-binary case**
      (`paths.namd3gpu == paths.namd3`), which the auto path could never elect — so a
      GPU-capable host with one binary serving both modes can now elect GPU via config instead
      of the fragile distinct-path-string workaround; `cpu` forces CPU. The decision is factored
      into a pure, unit-tested `Config._resolve_namd_type`. The `run --gpu` CLI flag remains as
      an equivalent per-invocation override. (Unreleased.)
- [x] **Offload all mmCIF parsing to pidibble.** mmCIF is now parsed by pidibble **≥1.7.1**
      (`PDBParser(input_format='mmCIF')`, which normalizes mmCIF into the PDB-record namespace);
      pestifer's raw-CIF layer (`util/cifutil.py` `CIFdict`/`CIFload`) is deleted and the direct
      `mmcif`/`mmcif-pdbx` dependencies dropped. Label-primary chain identity preserved; the only
      parsed-object behavior change is symmetry operators normalizing `1_555`→`1555` (aligns CIF with
      the PDB path). Two downstream fixes were needed (found running the examples): (a) the psfgen task
      hands the source `.cif` to VMD via `mol new`, whose pdbx plugin mis-parses a
      `pdbx_audit_revision_item` loop — `CIFload` used to strip it as a side effect, restored as the
      dependency-free `util.util.strip_cif_category`; and (b) mmCIF biological-assembly chain lists
      must be in the *label* namespace (pidibble's REMARK.350 `header` is author-mapped and lossy), so
      pidibble 1.7.1 adds `header_label` (the raw label `asym_id_list`) and `Transform` prefers it.
      Validated by a golden AsymmetricUnit oracle on 4zmj/6pti (diffs to only the sym delta), the full
      unit suite, and an all-examples runthrough (the CIF/glycan/assembly builds — 7–12, 14, 15, 18,
      21 — all pass). Design doc: `docs/design/mmcif-offload.md`.
- [x] **Ledger for `modify-package`.** Append-only record at
      `pestifer/resources/modifications.jsonl`; every mutating command records an entry
      (committed alongside the change). `ledger show` lists them; `ledger revert <id>`
      git-reverts a recorded modification on a fresh branch and curates the ledger.
      (v3.2.0.)
- [ ] **Restart / resume an interrupted build.** A long multi-task build that dies partway —
      Ctrl-C, a crash, a killed job, a wall-clock timeout — currently has to start over from
      scratch. Add a restart path that resumes from the last completed task instead of re-running
      the whole pipeline. Much of the substrate exists: each task already emits a state fileset
      (`psf`/`pdb`/`coor`/`xsc`/`vel`) as an artifact, the `continuation` task already knows how to
      begin from such a fileset, and the SIGINT/SIGTERM handler now tears children down cleanly on
      interrupt (v3.11.1) — the partial files it warns about are exactly what a restart would key
      off. Decisions to settle:
    - **Checkpoint granularity.** Per-task (resume at the first task whose outputs are missing or
      stale) vs. within-task (finer, but not every task is idempotent) — start per-task.
    - **Completion detection.** How to know a task finished cleanly and its outputs are current —
      e.g. a per-task completion manifest keyed by a hash of the task's resolved specs, so a config
      edit upstream of the resume point correctly invalidates and re-runs the tail.
    - **Invocation.** `pestifer run --restart` (auto-detect the resume point) vs. an explicit
      `--from <task>`, and how it composes with the existing `continuation` task.
    - **Determinism.** Lean on the seeded-RNG work (declash, CCD closure) so a resumed build is
      reproducible — ideally identical to an uninterrupted run from the resume point onward.
- [~] **Optionally-interactive `new-system` for sequence modifications.** `new-system` generated a
      build config from a PDB/UniProt ID off a fixed template with no look at the structure. Now it
      inspects the source for **missing residues** (chain gaps / `REMARK 465`) and **engineered
      mutations / tags** (`SEQADV`) and annotates the emitted YAML with the corresponding mods.
    - [x] **P1 — discovery + annotate (Unreleased).** `new-system <id> --inspect` fetches and parses
          the structure (pidibble, no `Molecule`) and surfaces: **biological assemblies** (index,
          copy count, chains — pick with `biological_assembly: N`, `0` = A.U.); **chain identities**
          (segtype + size, omit via `exclude:`); **missing residues** grouped into runs and
          classified interior-gap vs N-/C-terminal tail; and `SEQADV` engineered mutations/conflicts
          and cloning/expression-tag runs. Interior gaps get an active `ligate` task; every other
          finding is a **commented** stub under the `psfgen` task with correct nesting and
          ready-to-uncomment shortcodes (`biological_assembly`, `exclude`, `terminal_tails`,
          `mutations` to revert, `deletions` to excise). Core in `core/system_inspector.py`;
          validated on 4zmj.
    - [x] **P2 — interactive walkthrough (Unreleased).** `new-system <id> --interactive` prompts per
          finding (which assembly / A.U.? omit this chain? build this tail? revert this mutation?
          excise this tag? add a `ligate` task?) and injects the choices as *active* config — a
          `source:` block (`biological_assembly` + `exclude`), a `sequence:` block, a `mods:` block,
          and a `ligate` task, ready to run. Prompt logic factored (`interactive_select` +
          `make_prompter`) so it is testable without a TTY; blank input takes the shown default.
    - [x] **Chain identity + interior-loop options (Unreleased).** Chain identities carry the
          molecule name from the header (`COMPND` / mmCIF entity), and interior missing loops are
          presented in the interactive walkthrough with a "build in full vs. short GGG stub"
          (`substitutions`) choice.
    - [ ] **P3 — canonical-sequence diff.** Align the modeled sequence against the canonical UniProt
          reference to surface substitutions not recorded in `SEQADV` (no reuseable primitive exists
          today beyond the AlphaFold *structure* fetch — needs a UniProt FASTA fetch + alignment).

- [ ] **Leave a capped chain break in place of an unbuilt interior loop.** Today an interior
      missing loop is *always* built (then closed by `ligate`); the only alternative is a shorter
      built stub (`substitutions`). A user may instead want the loop *not* built and the two
      flanking runs left as **separate, properly-terminated segments** (a genuine chain break, each
      end capped with `CTER`/`NTER`). Validated that the obvious config-only route does **not** work:
      `exclude`/`deletions` drop the loop residues but *fuse* the anchors into one segment (an
      unphysical long bond), and a `cleave` task does not split at the resulting numbering gap
      (`ResID` is also not int-comparable, so a resid-range `exclude` fails outright). This needs a
      dedicated build-engine capability — split a chain at an unbuilt interior gap into two segments
      and apply default terminal patches — after which `new-system --interactive` can offer it as
      the third interior-loop option (build / stub / cap-and-break).

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

## Loop modeling

- [x] **Physics-based loop modeling for missing internal segments.** Replaced the
      `guesscoord` + steered-MD (`ligate`) closure — which yielded topologically-correct
      but stretched, unrealistic loops — with an offline, in-house **sample → CCD-close →
      score → minimize** pipeline. Fully self-contained: no structure predictor, no
      network. Shipped in v3.9.0 as the default `ligate method: ccd`; the `connect`
      peptide-bond patch is applied after closure and relaxation is deferred to a
      downstream `minimize`. `steer` demoted to an opt-in fallback (`method: steer`).
      A mid-project **scope correction** (see `docs/design/loop-modeling.md`) concluded
      the real workload is floppy, solvent-exposed surface loops that need only "no new
      clashes" rather than native reconstruction — so P2/P3 below were reframed as
      unnecessary for that goal. Full plan and history in `docs/design/loop-modeling.md`.
  - [x] **P1 — CCD closer replaces steering.** Analytic Ramachandran basins for the
        initial φ/ψ, cyclic-coordinate-descent closure onto the downstream anchor,
        `connect`, minimize downstream. Deterministic (seeded). `steer` → opt-in
        fallback. Shipped as-default in v3.9.0, and already does per-loop clash-filtered
        ensembles + iterative declash refinement (beyond the minimal P1). Engine in
        `psfutil/loop_ccd.py`; delete-and-rebuild BPTI benchmarks in the unit suite.
  - [ ] ~~**P2 — derived coil torsion library + ensemble ranking.**~~ **Dropped** per the
        scope correction: the analytic `RAMACHANDRAN_BASINS` in P1 suffice for surface
        loops, so a pestifer-generated coil φ/ψ library was not needed. Revisit only if a
        genuinely buried long loop demands native-quality reconstruction.
  - [ ] ~~**P3 — KIC + neighbor-dependent library + optional restrained-MD refine.**~~
        **Retired** — KIC was prototyped then dropped in the scope correction; the `steer`
        fallback is deliberately kept (not removed) as the opt-in escape hatch. Revisit
        only alongside P2.
  - [x] **Follow-up: missing terminal tails.** No downstream anchor, so no closure — a built
        terminal tail is now modeled by *sample + declash* (see the regularization item below):
        the psfgen free-tail modeler replaces guesscoord's extended arm with a Ramachandran-seeded,
        clash-filtered backbone, rigidly anchored at its resolved junction. (Unreleased.)

- [x] **Regularize how missing residues are modeled across terminal and internal segments.**
      Missing residues were handled by two divergent paths: interior gaps closed through the
      `ligate`/CCD machinery (build the run, close onto the downstream anchor, `connect`, minimize),
      while missing N-/C-terminal tails, when built at all, kept guesscoord's arbitrary extended arm
      with none of that quality machinery. Both halves now shipped (Unreleased):
    - **Modeling parity.** The interior closer's conformation engine (Ramachandran seed +
      clash-filtered ensemble + iterative declash) is shared with a terminal **free-tail modeler**
      (`psfutil/tail_model.py`) that dispatches by position — anchored interior gaps → CCD closure
      (`ligate`); terminal tails → sample + rigid junction placement + declash (psfgen declash step,
      `loops.declash.model_tails`) — sharing the seeded RNG, NeRF junction geometry, clash report,
      and the downstream `minimize`.
    - **Declaration parity.** The three flat opt-in keys are consolidated into one grouped
      `terminal_tails: {n, c, all}` block (`normalize_terminal_tails`), with the legacy keys folded
      in behind a one-time deprecation warning. Interior gaps stay always-built — the intended
      asymmetry (terminal tails are often disordered ends / expression tags to drop), so "declare
      the same way" means one grouped terminal surface, not forcing interior/terminal identical.
    - Composes with the optionally-interactive `new-system` item (which can emit detected missing
      residues straight into a `terminal_tails` block). Possible follow-up: ring-piercing hygiene
      for tails (interior loops get it via the graft-time scan).

## Ideas / unsorted

- [x] **Persist chain IDs across PSF→PDB regeneration.** A PSF has no chain column, so
      every `coor`→`pdb` regeneration re-derived each atom's chain from its segid's leading
      character — why merged copies collapsed onto one chain (v3.7.0 worked around it by
      forcing single-character segids in `merge`). `coor_to_pdb` now accepts a reference PDB
      and restores the chainID column from it by `(segid, resid)`
      (`restore_chain_ids_from_reference`); the MD task passes the pre-MD state PDB, so
      arbitrary chain IDs survive every MD regeneration of a build. The `merge`
      single-character-segid workaround is left in place (still the mechanism that *assigns*
      the distinct chains this restore then preserves); relaxing it is a possible follow-up.
      (Unreleased.)
- [ ] _(add items here)_
