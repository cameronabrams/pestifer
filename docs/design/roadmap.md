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

- [~] **Fluid-bilayer-like lipid conformers for grid packing.** The grid packer stamped one extended
      vacuum conformer onto every lipid, so gridded membranes started **over-packed / gel-trapped**
      (DMPC settles to APL ~50 vs. fluid ~60) and the slow post-embed area drift was the acyl chains
      slowly acquiring configurational entropy — a chicken-and-egg the packer couldn't escape. Fix:
      generate *fluid-like* conformers (moderately splayed, gauche-rich tails) with cylinder-confined
      **athermal Monte Carlo** and draw them per-lipid. Full design + results + the `autobonds off`
      dead-end in `docs/design/lipid-conformer-generation.md`. (All Unreleased.)
  - [x] **Lever 1 — per-lipid conformer draw.** `write_grid_pdb` caches the whole cached ensemble per
        species and draws a conformer per lipid (and per chamber solvent) from the packer RNG; explicit
        `conf` still pins. Fixed a latent `conformers` vs `lipid_conformers` key mismatch that had been
        silently ignoring the pin knob; schema default `'0'`→`''`. (`f8772a04`.)
  - [x] **Lever 2a — cylinder-confined athermal-MC conformer generator (opt-in).** New
        `charmmff/athermal_mc.py`: dihedral-pivot MC (bonds/angles preserved exactly; only torsions
        move), hard-sphere excluded volume (sigma-scaled `Rmin/2`, exclude ≤1-4 so gauche is allowed),
        monotone-inward cylinder confinement about the membrane normal through the tail-bundle centroid.
        Cylinder cross-section = `cylinder_inflation × APL` (a *packing* knob — the equilibrated APL is
        set by NPgT, not the conformers). Wired as `sampler='mc'` through
        `do_psfgen`/`do_resi`/`ensure_lipid_conformer`; **default stays `'md'`** until validated across
        lipid types. (`92d4953d`, `44fc6029`, `bd68dffa`, `3c242eff`.)
  - [x] **Sterols → single conformer.** Cholesterol-substream resis are forced to `sampler='single'`
        (init+minimize+orient, one rigid conformer); provenance records it. (`ccdc662e`.)
  - [x] **Packer robustness for fat conformers.** (1) A protein-free pristine bilayer was unbuildable
        from config — the schema always fills `embed` so the contract saw `embedding=True`; now detects
        a *meaningful* embed (`58b75144`). (2) Fluid conformers pack dense enough to (a) make VMD
        distance-bond near-coincident lipid atoms → **merge residues** → psfgen guesses → NaN, and (b)
        **segfault** the minimize on overlap; both fixed and *hardened* by the packer re-spin guard —
        best-of-40 attempts, escalating jitter + ensemble redraw, loud warnings, and a hard-abort below
        a true-coincidence floor (`c56863db`, `80326565`). (`autobonds off` was tried and reverted
        `9a316ed5` — the VMD split step needs bond perception; it exploded psfgen 23.6k→2.8M atoms.)
  - [x] **Acceptance test met (DMPC).** Clean control, identical protocol: old rods gel-trap at APL
        50.3; MC conformers reach fluid APL ~67 with 0.5% area drift — the over-packing bug is fixed.
  - [~] **Remaining:** (i) **ex16 payoff — DONE.** ex16 rebuilt end-to-end under the new builder
        (fluid MC conformers + solvate-based hydration): quilt converged step 180090, post-embed
        `membrane_equilibrate` converged step 288090, final DMPC density 1.0201 g/cc at APL ~61.2 Å²
        (spot-on for DMPC/310 K), all tasks result 0. The solvate hydration removed the water-void
        distortion and the fluid start cut the post-embed cost. (ii) **re-tune `membrane_equilibrate`
        area defaults** for the fluid start; (iii) **multi-lipid smoke test** (POPC/POPE/POPS/POPG, a
        sphingolipid, a mixture); then (iv) **flip `sampler='mc'` to default.** Deferred: absolute-APL
        vs literature (force-field / conditions, not conformer generation); Lever 2b (whole-grid MC
        pre-relax — likely unnecessary).

- [~] **Per-leaflet phase (Lo/Ld) as a build-time input.** Building ex17's 43–47%-cholesterol leaflets
      showed that *fluid* (Ld) conformers cannot calibrate a **raft** leaflet: a cholesterol-rich
      saturated leaflet undergoes a slow **liquid-ordered (Lo) phase transition** that
      `membrane_equilibrate` can't converge — the area drifts 48→46.7→44 Å² and consumes the whole 800k
      -step budget without plateauing (it's a phase change, not equilibration, so there's nothing to
      stop *at*). Since MD can't cross Ld↔Lo on build timescales, the phase you *build* is the phase you
      *keep*, so declare it up front and pack with the matching conformer ensemble instead of relaxing
      into it. Full spec in `docs/design/lipid-conformer-generation.md` (§ "Phase (Lo vs Ld) as a
      build-time input"). This is the root-cause fix for the calibration-reliability item below.
  - [x] **Ordering mechanism = a tunable trans-bias, NOT the cylinder.** Empirically disproved the
        first plan (tighten the cylinder to order the chains): an athermal (repulsive-only) MC stays
        fluid however tight the cylinder (chain order saturates ~0.15–0.23; the cylinder caps
        footprint/APL, not order — nothing drives the chains trans). Ordering is a separate *optional*
        trans-favoring torsional Metropolis term in `run_mc` (`torsion_bias`; potential `(1+cos φ)/2`,
        O(1)/move since a pivot moves exactly one dihedral; `bias=0` reproduces the athermal path
        bit-for-bit). **Ld = bias 0** (the athermal fluid ensemble *is* Ld); **Lo = bias auto-tuned**
        (`_bisect_to_order`) to a chain order parameter (Scd) target measured by
        `chain_order_parameter`/`ensemble_chain_order` (−mean ½(3cos²θ−1) over tail C–H vs +z). Unit-
        tested. (`3b641833`.)
  - [x] **Per-leaflet integration.** schema `composition.upper/lower_leaflet_phase` (choices `[Ld, Lo]`,
        default `Ld` = prior behavior); the packer routes conformer lookup through a per-species
        `conf_key` (`PSM__Lo`) while the *placed* resname stays the real CHARMM resname (it comes from
        the conformer file, not the key); phase-aware autocache entries cache as `<resname>__<phase>`
        without colliding with the shipped fluid set; a phase/sterol composition-mismatch warning.
        (`fa2890ae`; salt-ion `conf_key` fix `c912994f`.)
  - [x] **Validated — the raft calibration now converges.** ex17's upper leaflet (PSM 36 / POPC 17 /
        CHL1 47) declared `Lo`: patch-A calibration **plateaus at APL ~46.8 Å² (converged step ~258k,
        area drift +0.002)** — versus the old Ld run's non-converging 48→44 drift that ate all 800k
        steps. Packing in-phase removes the Lo transition from the critical path.
  - [x] **The unsaturated-Lo ceiling is self-correcting (not a bug).** The trans-bias can't order
        *unsaturated* chains — the *cis* double bond is a rigid kink, not a rotatable dihedral — so
        POPE/SOPS/SOPE clamp at the bias ceiling (order ~0.22 ≈ fluid) while saturated PSM/POPC reach
        0.31–0.33. That is *physically correct*: saturated+cholesterol makes a liquid-ordered raft,
        unsaturated+cholesterol stays fluid. So declaring an unsaturated leaflet `Lo` harmlessly yields
        the fluid ensemble it *should* have. Confirmed on ex17: patch-B (SOPE/SOPS/POPE + 43% CHL1)
        **converged cleanly at a fluid APL ~56.0** (no ordering drift), while patch-A (saturated raft)
        converged dense at ~46.8 — a real asymmetric raft/fluid bilayer. Follow-up (efficiency, not
        correctness): declare such leaflets `Ld` to skip the wasted Lo generation, since Ld and Lo give
        an unsaturated leaflet the same fluid ensemble.
  - [ ] **Remaining:** calibrate the Lo Scd target against *packing* (does order 0.30 hold the dense
        basin, or is the per-lipid conformer too weak a lever?); skip redundant `<sterol>__Lo` entries
        (sterols are phase-independent, forced `single`); finish the ex17 both-`Lo` acceptance build
        (patch-B → quilt → embed → equilibrate) before committing the example.

- [ ] **Decide whether `max_steps` is a hard ceiling or a soft budget.** A chunked
      equilibration runs whole chunks, so the last one crosses the ceiling rather than stopping at
      it: a run with `max_steps: 6000` was observed to perform 8000 steps and then report
      `CEILING: reached max_steps (6000)`. Both numbers are true -- the run *did* exhaust its budget
      and it *did* run 8000 steps -- but a reader comparing them will assume one is wrong. Either
      quantize the final chunk down to the remaining budget, or say "budget" rather than
      "max_steps" in the message and the schema. (Observed while writing the loop tests, 2026-08-19.)

- [ ] **Compute and plot chain order parameters (Scd) alongside density / box size.** The phase work
      packs and calibrates leaflets by chain order, but that order is only ever measured on the
      *single-molecule conformer ensemble* at generation time -- never on the assembled, equilibrating
      membrane. Make Scd a first-class diagnostic, plotted the way density and box dimensions already
      are: a **mean-Scd time series** (convergence) and the standard **per-carbon `|S_CD|` vs carbon-
      number profile**; add **density(z) profiles** where applicable. The numeric core already exists
      (`athermal_mc.chain_order_parameter` / `tail_carbon_indices` / `ensemble_chain_order`, built for
      the phase work); the remaining work is coordinate access + plotting.
  - [ ] **P1 -- Scd in the `membrane_equilibrate` convergence plot (moderate; no new dependency).** That
        task already tracks density + lateral area per chunk and emits `-membrane.png`, and each chunk
        leaves a `.coor`/`.pdb`. Read the chunk coords + PSF bonds (`PSFContents`), compute Scd with the
        existing helper (extended to per-carbon, ~20 lines), and add an Scd-vs-time panel plus a final
        per-carbon `|S_CD|` profile -- showing the raft leaflet *ordering* alongside density/area, the
        exact quantity this arc targets.
  - [ ] **P2 -- Scd + density(z) profiles on the production trajectory (harder).** mdplot reads only
        scalars from the NAMD log (density, cell) and NAMD's own pressure profile; it never touches
        coordinates. Scd and density(z) both need the DCD, which pestifer currently handles only via VMD
        + `catdcd` (no Python trajectory reader). Decide: a VMD script computing/dumping per-frame, vs.
        adding a lightweight trajectory reader (MDAnalysis would make this -- and any future coordinate
        analysis -- much easier, at one dependency). Density(z) profiles come nearly for free once
        trajectory reading exists.

- [x] **Reconcile the two membrane area-convergence criteria on the calibration path.** *Done.* When an
      asymmetric build's symmetric **calibration patch** relaxes via a self-terminating
      `membrane_equilibrate` (migrated from the old fixed NPgT ladder), two different area-convergence
      tests now both apply: `membrane_equilibrate`'s own autocorrelation-corrected gate (default
      `area_drift_tol = 0.012`, a deliberately loose soft-mode default) and `make_membrane_system`'s
      post-hoc `_area_convergence` final-half-drift reliability check (`_AREA_CONVERGENCE_TOL = 0.005`).
      They can disagree — a patch converges under the loose gate while still creeping at ~0.7% and the
      stricter check then warns the calibrated preferred APL is unreliable (observed on ex17's
      47%-cholesterol upper leaflet: converged APL ~45.4 with −0.67% residual drift). Worked around per
      example by tightening the patch stage's `area_drift_tol` to 0.005, but the redundancy is a design
      smell: the `_area_convergence` re-check was built for the old *non*-self-terminating fixed-ladder
      patches. Make `make_membrane_system` **trust a converged `membrane_equilibrate`'s own verdict**
      (skip/relax the cruder re-check when the last patch stage self-terminated), or unify the two on a
      single autocorrelation-corrected criterion, so calibration precision is set in one place.
      **Resolved by the first option.** `membrane_equilibrate` gained an opt-in cumulative
      `area_plateau_tol` gate measuring the *same* trailing quarter-over-quarter drift that
      `_area_convergence` uses, and `make_membrane_system._relaxed_area_drift` now trusts a
      plateau-gated stage that CONVERGED, returning its own plateau drift and falling back to the
      post-hoc check only for an `md`-terminated ladder or a stage that hit its ceiling. The two can
      no longer disagree, and the ex17 workaround (hand-tightening `area_drift_tol`) is retired in
      favor of `area_plateau_tol: 0.005` on the calibration protocol.

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
- [ ] **Auto-resolve sterol IC-donor dependencies in conformer generation.** `CHM1` (and possibly
      other sterols) has an intentionally empty IC table in the bundled
      `toppar_all36_lipid_cholesterol.str` (bond lengths/angles all `0.0000`): it is *designed* to
      inherit its geometry from `CHL1`, a normal CHARMM cross-residue IC pattern — not a typo. psfgen
      coordinate guessing therefore fails during on-the-fly generation only because pestifer doesn't know
      to supply the donor; it's currently worked around by hand with `--take-ic-from CHL1` (CHM1 is a
      `single`-conformer sterol, so no MC is involved). Ship a small **known IC-donor map** (or auto-
      detect an all-zero IC table and fall back to the parent-ring donor) so generation resolves these
      dependencies itself instead of failing and needing a manual retry.
- [ ] **Audit sterol RESI `BOND` records for mis-assignments.** Separate from the IC-donor issue: at
      least one sterol RESI block (recollection: CHM1 or another in the cholesterol substream) is
      believed to have **incorrectly assigned bonds** in the vendored CHARMM stream, which would give a
      wrong psfgen topology and a distorted built conformer. Audit the cholesterol-substream RESI blocks'
      connectivity against the intended sterol skeleton; if confirmed, ship a committed **diff patch**
      against the bundled stream. Fits the general "carry local patches over vendored CHARMM streams"
      need — see the CHARMM-release ingestion item below, which should reapply such patches on each
      release bump.

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

- [ ] **The per-user cache is keyed on version but not on which checkout wrote it.**
      `CacheableObject` names its files `<name>-v<major>.<minor>`, so two checkouts at the same
      version share one cache -- and the cached CHARMMFF content stores *absolute paths* to
      resource files. Building from a second checkout (a git worktree, a pip-installed copy
      alongside the source tree) therefore poisons the shared cache with paths into that copy, and
      builds from the original then fail with `FileNotFoundError` on a toppar file. `pestifer cache
      clear` fixes it, but nothing points there. Options: include a hash of the resource root in the
      cache key, or store resource paths relative to it. (Hit while running the test suite from a
      worktree, 2026-08-19.)

- [ ] **`deposit` subcommand: package a completed build for public archiving.** Point it at a
      finished production fileset and it assembles a deposit-ready bundle — inputs, topologies,
      parameters, initial configurations, and outputs — laid out for upload to Zenodo or an
      equivalent DOI-issuing repository. This maps directly onto what JCIM now requires of an MD
      paper: *"depositing your input files (including topologies, parameters, input files, and
      initial configurations) and output trajectory files in a public repository providing a
      DOI."* Points to settle:
    - **Trajectory size.** Full trajectories will usually be too large — 27 examples × 3 replicas
      is already substantial. The guidelines permit depositing *a selection of clustered
      representative structures* instead, so an option to cluster and deposit representatives
      (rather than everything) is probably the common path, with full-trajectory deposit as the
      opt-in.
    - **What it reuses.** A build already writes `run-record.json`, an environment record, and a
      citation list; the deposit should carry those rather than re-derive them, so the bundle
      documents its own provenance. It pairs naturally with `report-methods`: one drafts the
      Methods text, the other assembles the deposit, and between them a build becomes publishable
      without hand-assembly.
    - **Sequencing.** Same argument as environment recording: a deposit assembled by a tool that
      already knows how the build was produced is reproducible by construction, whereas one
      hand-assembled afterward is a reconstruction. Worth having *before* a deposit is made in
      anger.
    (Raised by Cameron while preparing the manuscript, 2026-08-18. Deliberately **not**
    implemented at the time: a new feature means a new release, and the manuscript pins v3.18.0
    as the version that produced its results.)

- [x] **GPU mode detection shouldn't hinge on path inequality.** The `namd.processor-type`
      schema entry (formerly deprecated-and-ignored) is now a functional `auto|cpu|gpu` control
      (default `auto`). `auto` keeps the historical path-inequality detection (+ the
      `_report_gpu_mode` status line); `gpu` **forces GPU mode even in the single-binary case**
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
- [ ] **Extract the NAMD analysis layer into its own package.** Reading and analysing what NAMD
      writes -- `.log` and `.xst` parsing, density and lateral-area computation, the
      autocorrelation-corrected convergence criterion, density profiles, restart generation -- is a
      general capability worth having outside pestifer, on the precedent of `pidibble` and
      `ycleptic`. Roughly 1,000 statements with light or no coupling back into pestifer.  *Writing*
      a NAMD config is not extractable: it needs the CHARMM force field and the PSF, so it is
      pestifer's own domain.  Note pestifer would **depend on** the result rather than shed it, as
      it does pidibble.  Full evidence, coupling graph, sequencing and the decision test in
      `docs/design/namd-layer-extraction.md`.  (Distinct from the resource-package item below.)

- [ ] **Split the stable CHARMM resources into a separate data package.** The bundled force
      field and PDB repository are *data*, versioned on CHARMM's release cadence, but they live in
      the code repo and ride every pestifer release. The costs, measured rather than supposed:
    - **PyPI quota — real but no longer pressing.** This was the original motivation, when the
      preflight in `scripts/release.sh` reported 5.91 GB of PyPI's 10 GB per-project limit
      consumed. Deleting the pre-3.7 releases (~100 MB each, from before `tests/**` was excluded
      from the sdist) brought that to **3.56 GB across 49 published versions**, and current
      releases are ~24 MB, so there is headroom for **~270 more** — not a wall. Two caveats keep
      it on the list rather than off: published files count forever unless manually deleted, and
      deleting them breaks reproducibility for anyone pinning an old version, so trading history
      for headroom has a cost. The force field is still the dominant share of each release
      (`pestifer/resources/charmmff` is 9.4 MB of a 25 MB package) and does not compress further —
      it is already tarballs.
    - **Git history — the strongest argument.** Force-field paths account for **212 MB of the
      818 MB packed history**, and **155 MB of that is just 88 blobs** — the
      `pdbrepository`/`toppar` tarballs, each of which
      is rewritten *whole* on every regeneration because a tarball has no delta. Rebuilding the
      lipid collection with a new sampler costs another full copy forever. Every clone pays it.
    - **Coupling.** A force-field fix and a code fix are the same event today: correcting a
      vendored stream (see the sterol `BOND` audit above) requires a pestifer release, and a
      pestifer patch release re-ships an unchanged 9.4 MB force field.

    The natural seam already exists. A release directory is self-contained and release-keyed
    (`charmmff/<key>/` = `toppar_c36_<key>.tgz` + `pdbrepository/` + `patches/`, keyed by
    `charmmff_version_key`), the config already selects a release (`charmmff.version`), and the
    autocache already treats `~/.pestifer/pdbrepository/<release>/` as a first-class *external*
    source of the same content. So pestifer already knows how to consume this data from
    somewhere other than its own package directory.

    Decisions to settle:
    - **Distribution mechanism.** A companion PyPI package (`pestifer-charmmff`) pinned as a
      normal dependency — simplest, keeps `pip install pestifer` a one-liner, but each release
      still consumes quota, just in a *separate* project with its own 10 GB. Versus fetch-on-first-
      use into `~/.pestifer/` from a data repo or Zenodo record — unbounded, content-addressable,
      citable (a paper can cite the exact FF release DOI), but makes the first run network-
      dependent and needs an offline/air-gapped story for HPC. The autocache's existing
      lock + atomic-publish machinery is the substrate for the second.
    - **Granularity.** One package per CHARMM release (`pestifer-charmmff-feb26`), letting a user
      hold several and switch with `charmmff.version`, versus one package tracking the latest.
      Per-release matches how the code already keys things and makes reproducing an old build a
      matter of installing an old data package.
    - **Where local patches live.** `patches/` and `custom/` are pestifer's own corrections and
      additions on top of vendored CHARMM, not upstream content. They probably stay with the
      code (they change on pestifer's cadence, and the CHARMM-release ingestion item below is
      supposed to *reapply* them across a release bump) — but then the data package is pure
      upstream, and the split has to survive `modify-package charmmff regenerate-*` and the
      `modifications.jsonl` ledger, which currently assume one writable tree.
    - **Generated vs. vendored.** `pdbrepository` is *derived* (conformers pestifer generated),
      while `toppar` is vendored. They have different provenance and different rebuild triggers,
      so they may deserve separate packages rather than one bundle.

    Composes with **"Automatic update to a new CHARMM force-field release"** below: that item
    automates *producing* a release directory, and this one decides where the product is
    published. Doing this one first would make that one's output a data-package release rather
    than a commit to the code repo.

- [x] **Ledger for `modify-package`.** Append-only record at
      `pestifer/resources/modifications.jsonl`; every mutating command records an entry
      (committed alongside the change). `ledger show` lists them; `ledger revert <id>`
      git-reverts a recorded modification on a fresh branch and curates the ledger.
      (v3.2.0.)
- [x] **Restart / resume an interrupted build.** *Done (v3.15.0).* A long multi-task build that dies partway —
      Ctrl-C, a crash, a killed job, a wall-clock timeout — currently has to start over from
      scratch. Add a restart path that resumes from the last completed task instead of re-running
      the whole pipeline. Much of the substrate exists: each task already emits a state fileset
      (`psf`/`pdb`/`coor`/`xsc`/`vel`) as an artifact, the `continuation` task already knows how to
      begin from such a fileset, and the SIGINT/SIGTERM handler now tears children down cleanly on
      interrupt (v3.11.1) — the partial files it warns about are exactly what a restart would key
      off. Every decision below was settled as proposed; see `docs/design/build-restart.md`.
    - **Checkpoint granularity** — per-task, as suggested.
    - **Completion detection** — a run manifest (`.pestifer-manifest.json`) records each task's
      spec hash, taken *before* execution since some tasks mutate their own specs, so editing the
      config invalidates from the first changed task onward.
    - **Invocation** — `--restart` auto-detects the resume point, `--from <task>` overrides it,
      `--fresh` ignores an existing manifest. A resume point must be a clean `STATE` hand-off;
      pestifer refuses to resume into a task whose predecessor produced no state.
    - **Determinism** — unchanged from the plan; rests on the existing seeded-RNG work.
      Documented for users under "Resuming an interrupted build" in the `build` subcommand page.
- [~] **`density_equilibrate` task — convergence-based post-solvation equilibration.** Replace the
      hand-written ladder of progressively longer NPT runs (`nsteps: 200 → 400 → … → 13200`) with a
      single task that runs NPT in chunks and stops itself when the **box density has converged** or a
      `max_steps` ceiling is hit. De-magic-numbers the workflow and adapts to system size (small boxes
      stop early, large ones run longer) instead of a fixed schedule. Reuses the existing MD
      continuation substrate (`firsttimestep` + state fileset) and reads density from the NAMD `.xst`
      cell volume + PSF mass; composes with the restart/resume item above. Seeded → reproducible stop
      on a fixed machine (not bit-identical cross-platform, by nature of wrapping NAMD). Full plan in
      `docs/design/density-equilibrate.md`.
    - [x] **P1 — the task (Unreleased).** `density_equilibrate` (subclasses `MDTask`) with the
          stability-bounded adaptive chunk sizing, the precision-gated fractional-drift criterion with
          `n_consecutive` hysteresis, `max_steps` guard and NaN abort, and a convergence-report +
          density-vs-time-plot artifact. **Reactive stability net:** each `namdrun` is wrapped in a
          patch-grid crash-catch — detect NAMD's "cell too small for patch grid" abort, roll back to
          the previous chunk's registered state (free — a crashed run registers none) and retry
          shorter with AIMD re-growth (halve on crash, grow ≤`chunk_growth`× per success), so
          correctness doesn't depend on modeling NAMD's exact threshold and CPU runs (smaller effective
          `margin`) self-tune instead of re-crashing every restart. Convergence + crash-detection +
          patch-grid parsing all live in a NAMD-free, unit-tested `util/density_convergence.py`. Wired
          into the `new-system` template and interactive pipeline (replacing the NPT ladder there).
          **Validated on real NAMD (BPTI, small-box stress case):** running the task end-to-end found
          three things unit tests could not — a non-`stepsPerCycle` chunk length (NAMD aborts; fixed
          with whole-cycle quantization), a full-history convergence window that never converges on a
          plateaued box (fixed with a trailing `window_frac`), and a `drift_tol` gate too tight for a
          15k-atom box (relaxed 0.1%→0.2%, calibrated over three independent trajectories via offline
          replay of the real monitor). The shipped defaults self-terminate at the plateau density (live
          run stopped at step 71530, before the ceiling). **Still to do (nice-to-have):** confirm on a
          large solvated box (expected to converge *earlier* — `SEM/mean`∝1/√N) and a membrane (keeps
          the NPgT protocol; verify the orthorhombic-volume density is sane).
    - [x] **P2 — migrate the bundled examples (Unreleased).** All 25 soluble examples replace their
          NPT ladder + `mdplot` with a single `density_equilibrate`; the 2 membrane examples (16, 17)
          stay on `NPgT` (out of scope). Comment-preserving text-level migration; all 27 configs still
          schema-validate. Validated end-to-end on **GPU** by building the insulin hexamer (example 13,
          ~2× BPTI, GPU-resident NAMD `margin: 4`): clean build, no patch-grid crashes, converged at
          step 50640 (< the 100k ceiling) at the box plateau — the criterion correctly waited out a
          higher/slower densification (→1.06 g/cc) without a premature stop. A truly large box
          (GroEL/HIV-Env, 100k+ atoms) is the remaining nice-to-have confirmation (expected to converge
          *earlier*, `SEM/mean`∝1/√N) but needs more GPU memory than the local 4 GB card.
          **Full 27-example CPU reference sweep: 27/27 pass, 0 failures** (large boxes up to ~1M atoms
          confirm the 1/√N trend on CPU; reactive net never needed). See P2.5 below for the one finding.
    - [x] **P2.5 — honest, size-aware convergence (Unreleased; supersedes the size-gate prototype).**
          The over-conservatism turned out not to need an ad-hoc size gate: a fine-sampling (`xstfreq=1`)
          experiment measured the density autocorrelation time (τ_int ≈ 559 steps — the block-means SEM
          under-reported uncertainty because `.xst` frames are correlated), so the precision gate now
          uses an **autocorrelation-corrected SEM** (`σ/√N_eff`, `N_eff = N/τ_int`) with a reliability
          guard. This is size-aware *for free* (`σ/mean ∝ 1/√N_atoms`, τ ~size-independent → big boxes
          converge in the fewest chunks; small boxes are fluctuation-limited), so the prototyped
          `(N_ref/N)^0.5` gate was dropped as redundant. The drift test became an **equivalence/upper-
          confidence-bound** test (fixing the plateau flicker without the √12 significance floor). A
          second full 27-example sweep validated it: **24/25 converge + 2 membranes = 27/27, 0 errors**;
          organic solvents improved (acetonitrile 49,850 vs. old 94,540; DMSO 53,740 vs. 64,470).
          **Residual:** only **acetone** still hits a *benign* `max_steps` ceiling (precision met, drift
          4e-4, density flat) — the one genuinely high-noise small-organic box; a **per-solvent
          `drift_tol` override** is the remaining low-risk follow-up. See `density-equilibrate.md`
          ("Honest autocorrelation-corrected criterion").
    - [x] **P3 — membrane-aware equilibration (density + lateral area) (Unreleased).** A
          self-terminating NPgT sibling of `density_equilibrate` that stops when BOTH box density and
          membrane lateral area (`A = a_x·b_y`) have converged (`V = A·c_z` pins the anisotropic cell).
          Built on the observable-agnostic engine (`xst_cell_areas` + `JointConvergence`) with a
          base-class refactor into shared `ChunkedEquilibrateTask` / `DensityEquilibrateTask` (NPT) /
          `MembraneEquilibrateTask` (NPgT, [density, area]). Shipped M1 (engine primitives, `a6621ae6`)
          → M2 (task + measured area autocorrelation, `31f46fc2`) → M3–M5 (`feda3740`: wired into
          `make_membrane_system` pre-embed quilt + examples 16/17, per-leaflet protein-corrected APL,
          soft-mode area criterion with `area_min_steps` floor, the **100/50 NPgT barostat**, and M5
          whole-box embed solvation). Validated end-to-end on ex16. Design doc:
          `docs/design/membrane-equilibrate.md`.
        - **Residual cleanup:** the old `_autostage_protocol` fixed-NPgT tail and the crude
          `_area_convergence` helper are still present/called in `make_membrane_system.py`
          (lines ~490/637/885) alongside the new task; retire them.
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
          presented in the interactive walkthrough with a "build in full vs. short stub"
          (`substitutions`, user-chosen sequence) choice.
    - [x] **Full pipeline builder + launch (Unreleased).** `--interactive` also walks the post-psfgen
          pipeline (vacuum min, vacuum MD, solvate, solvated min, NVT/NPT equilibration, production,
          mdplot, `terminate`) so the whole YAML is assembled interactively, then offers to launch
          the build immediately (`interactive_pipeline`).
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

- [ ] **`--kick-ass` loses the banner's provenance lines.** The hidden `--kick-ass` flag swaps
      the plain banner for the block-art logo, but `_enhanced_banner_message` carries no
      `pestifer_version` or `charmmff_version` placeholders, so a run logged with it drops the
      two facts the ordinary banner records: which pestifer built the system, and which CHARMM
      force-field release it drew on. That makes the flag a trade of provenance for style, when
      it should be a trade of nothing. Setting the version and force-field lines beneath the
      logo would make the two banners interchangeable.

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

## Testing / CI

- [ ] **Give the integration tests somewhere to actually run.** They build real systems and are the
      only tests that can catch a wrong *structure*, but GitHub-hosted runners have no VMD or NAMD,
      so adding them to `.github/workflows/tests.yaml` would make them **skip** and report green --
      worse than not running them, because it looks like coverage. They currently run in
      `scripts/release.sh` (gated, ~2 min) and nowhere else, which means a regression is caught at
      release time rather than at the commit that caused it. A self-hosted runner on a machine with
      the toolchain is the real fix. See `docs/source/contributing.rst`.

- [x] **Coverage can be measured, and two blind spots are gone.** `pytest-cov` was not installed,
      so the distribution table could not be produced at all. With it: **72% overall** (23,917
      statements). The shape is healthy -- the data model is well covered (`objs` 86%, `molecule`
      85%, `psfutil` 83%) and the gaps are in orchestration (`tasks` 54%, `subcommands` 47%),
      most of which only executes with VMD and NAMD present.

      Two modules were brought up from near nothing: `util/namdrestart.py` 13% → 99% (which
      surfaced a real bug -- see the CHANGELOG) and `sphinxext` 38% → 95%, with `task_table` and
      `tclscript` both at 100%. The `sphinxext` figure had been *fictional*: its `conftest.py`
      skips the directory when sphinx is absent and sphinx was not a test dependency, so 15
      existing tests ran nowhere. That is the same failure mode as the integration-test item
      above, and it is worth watching for elsewhere -- **a conditional skip that is always taken
      reports green and looks like coverage.** (Unreleased.)

- [x] **The August coverage gaps are closed.** Every module flagged then has been brought up
      without a toolchain:

      | Module | Before | After |
      |---|---|---|
      | `tasks/mdplot.py` | 42% | 90% |
      | `tasks/make_membrane_system.py` | 23% | 54% |
      | `tasks/psfgen.py` | 47% | 65% |
      | `charmmff/make_pdb_collection.py` | 14% | 43% |
      | `tasks/ringcheck.py` | 20% | 99% |
      | `util/namdrestart.py` | 13% | 99% |
      | `sphinxext` | 38% | 95% |

      The lesson, which held in every case: **"needs the toolchain" was about half true.** Grid
      placement, psfgen and NAMD genuinely cannot be stubbed, but what each task *decides* around
      those calls is ordinary logic -- and it is the half where errors are silent, because a
      toolchain failure announces itself while a wrong decision just produces a structure that
      builds and runs. Where a task really is bound to the toolchain (`do_psfgen`, tail modeling,
      the numpy declashers), that is now the honest remainder rather than an assumption.

      Two real defects came out of the exercise: the `make-namd-restart` basename bug, and the
      `sphinxext` tests that had never run anywhere. (Unreleased.)

      `psfgen` is now 47% → 65%; what is left there is genuine structural work (tail modeling, the
      numpy declashers over real PSFs, molecule ingestion) rather than decision logic.
      (Unreleased.)

      `tasks/make_membrane_system.py` is done as far as it usefully goes without the toolchain:
      23% → 54%. The remainder is genuinely VMD/NAMD-bound (grid placement, psfgen, the relaxation
      ladder, protein embedding). The lesson worth reusing is that a task being "toolchain-bound"
      is usually about half true -- what it *decides* before and after shelling out is ordinary
      logic, and it is the half where errors are silent. (Unreleased.)

      `tasks/mdplot.py` is done: 42% → 90%, which was the largest single gap. The approach worth
      reusing on the others is the third layer of it -- driving the task's `do()` through a real
      `Config`/`Controller` against a **trimmed real** NAMD log, which needs no VMD or NAMD and
      cost ~4 s. Anything that consumes NAMD output rather than producing it can be tested that
      way; the fixtures are ~30 KB each. (Unreleased.)

- [ ] **Test runs dirty the working tree.** *Half fixed.* Several files under
      `tests/unit/test_tasks/test_psfgen_*/` and `test_cleave/` are tracked *and* rewritten by the
      suite, so `git status` is dirty after every run -- noise at best, and at worst it trips
      `release.sh`'s clean-tree precondition.

      The **timestamp** cause is gone: the generated-script header now carries `pestifer <version>
      seed <seed>` instead of `Created <ctime>` (see *Watermarks on generated files* in
      `build-provenance.rst`), and `test_psfgen_preserve`'s `.tcl` and `test_cleave`'s
      `_minimal.prm` are now byte-stable across runs.

      What remains is **a different cause**, found while verifying the above:
      `test_psfgen_resegment/00-01-000_psfgen-build.tcl` still churns, because the generated script
      embeds process-global counters -- `set m13 ...` vs `set m8 ...`, `Transform 14 begins` vs
      `Transform 6 begins`. Those depend on how many molecules and transforms were created earlier
      *in the same interpreter*, so the file's content depends on which other tests ran first.
      Note this does **not** weaken the byte-reproducibility of a build: for a fixed task list the
      counter sequence is fixed. It is test-ordering sensitivity, and the fix is the same as
      before -- write to a temp directory, or stop tracking generated output and keep only the
      genuine input fixtures.

- [x] **Cover the chunked-equilibration machinery.** The convergence *criterion* was well tested
      (94%) but the loop around it was not (23%): chunk sizing and growth, the patch-grid
      crash-halve-retry path, budget accounting relative to an inherited `firsttimestep`, and the
      stop reasons that end up in `run-record.json`. Both `ChunkedEquilibrateTask` and the density
      and membrane hooks are now driven without NAMD, by scripting the verdicts and writing
      synthetic `.xst` series. `equilibrate_base` 23% → 80%, `density_equilibrate` 21% → 46%,
      `membrane_equilibrate` 42% → 73%. (Unreleased.)

- [x] **A shared structural battery for integration tests.** `tests/integration/helpers.py`
      supplies `assert_psf_sane`: PSF and PDB describing the same system, no duplicate/self/
      over-long bonds, nothing stranded at the origin, integral total charge, no atoms occupying
      the same space. Each check corresponds to a bug that has actually shipped here. New tests for
      loop closure, cleavage, chain identity and bilayer packing use it. (Unreleased.)

**Deliberately not planned:** image-comparison tests for the plotting code. Baseline PNGs are
hostage to matplotlib, FreeType and font versions, answer "something changed" rather than "what
changed", and bless whatever bug was present when the baseline was generated. `mdplot` is tested
instead through its pure data-preparation functions and by asserting on the *content* of the figure
object (lines drawn, labels attached, axis text) -- which is how the blank-figure defect was found.
Much of the remaining uncovered plotting code is styling, where the only honest assertion is "it
did not crash"; testing it would move the coverage number without improving the odds of catching a
real bug.

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
- [x] **Bring an incoming (foreign) PSF into the pipeline (build *onto* an existing topology).**
      *Complete (Unreleased): P1 + P2.0–2.5 + P3 all shipped and pushed.*
      Take a **pre-built system** (CHARMM-GUI, a prior pestifer build, another tool) into pestifer's
      downstream machinery *without rebuilding*, and later *extend* it (ligand, glycan, mutation,
      patch). **Architecture decided:** state enters through a task that *provides* the `STATE`
      currency — never a `psfgen` source key (psfgen consumes `SOURCE` and builds; it is not a state
      door). The contract-correct entry is **`continuation`**, which already reads an external
      psf+pdb/coor(+xsc/vel) and even resolves a foreign PSF's recorded stream files. Full design:
      `docs/design/incoming-psf.md`.
  - [x] **P1 — ingest + classify + force-field-consistency preflight (Unreleased).** A
        `verify_parameters` flag on `continuation` (default True; pestifer's own internal
        continuations set it False to keep the zero-parse fast path) parses the incoming PSF once to
        (a) report a per-segid composition inventory and **warn on resnames unknown to the segtype
        table**, and (b) **hard-error up front** if the build's CHARMM release cannot resolve every
        atom type and bond/angle/dihedral/improper term — naming the residue and term — instead of a
        cryptic mid-NAMD `DIDN'T FIND vdW PARAMETER` abort. New NAMD-free, unit-tested checker
        `charmmff/psf_param_check.py` (CHARMM-faithful wildcard matching; dedup by type-tuple → cheap
        on huge systems; **zero false positives on a real BPTI PSF vs. full feb26**). `PSFAtom`
        resname→segtype lookup made non-crashing (`.get`) so a foreign custom-ligand resname
        classifies blank rather than `KeyError`-ing the parse. Pass-through interop
        (`continuation` → `md`/`membrane_equilibrate`/`terminate`) works today.
  - [~] **P2 — additive edits.** A downstream `psfgen` consuming the `STATE` via a **readpsf-preserve**
        mode (`load_project`, `guesscoord=False`/`regenerate=False`), *not* segment-rebuild.
        - [x] **P2.0 — plumbing (Unreleased).** Pipeline-aware contract (`validate_pipeline` passes
              `available` to contracts that declare it; `psfgen` infers preserve-vs-build from what
              precedes it — STATE + no SOURCE → preserve, else build) + `PsfgenTask._psfgen_preserve()`
              readpsf pass-through; mods hard-error as "not yet supported". `continuation → psfgen`
              validates and reproduces the system atom-for-atom; build path unregressed.
        - [x] **P2.1 — patches (Unreleased).** `patches` emitted as psfgen `patch` commands on the
              loaded project + `regenerate`/`guesscoord`; supported-mod allowlist hard-errors on the
              rest. Verified: ASPP protonates ASP A:3 on an incoming BPTI PSF (+1 atom).
        - [x] **P2.2 — ssbonds (Unreleased).** `ssbonds` route to `patch DISU` on the loaded project;
              allowlist now `{patches, ssbonds}`. Verified on BPTI (routing/acceptance/clean run).
        - [x] **P2.3 — coord torsion rotations (Unreleased).** `irotations`/`crotations` apply via
              `coormods()`; the preserve path builds the base molecule lazily first
              (`ensure_base_molecule` → identity assembly) so the per-image rotation has its transforms.
              Verified: CHI1 on ASP A:3 of an incoming BPTI PSF pivots the sidechain (CB fixed, OD2
              swings). Establishes the lazy-molecule pattern P2.4 reuses.
        - [x] **P2.4 — links (Unreleased).** A `links` mod's patch is resolved from residue geometry
              (`assign_residues` → `set_patchname` → `ic_reference_closest`) using the lazy molecule,
              then emitted via `write_links` (identity assembly transform). Verified on a committed
              glycoprotein fixture (ASN61–glycan resolves to NGLA, emits `patch NGLA A:61 V:1304`).
        - [x] **P2.5 — grafts, additive only (Unreleased).** A graft *adds* a glycan into a fresh
              segment (`segment{}` + `coordpdb` + cross-link patch) via `_emit_grafts_preserve`. Only
              **additive** grafts fit: extending a *terminal* receiver. A graft onto an internal
              residue (has downstream) would *replace* it (remove residues = re-segmenting) and
              hard-errors as P3. Fixed an `atom.elem`-unset mismatch in the aligner for prebuilt bases.
              Validated on `cleave_inputs` + a committed `mannose_donor.pdb` (extend terminal V:1314;
              internal V:1301 → P3 error). Fetch-metadata mods (`biological_assembly`, `SEQADV`,
              `REMARK 465`, `terminal_tails`) and re-segmenting mods (mutations/deletions/insertions)
              still hard-error.
  - [x] **P3 — re-segmenting edits: mutations / deletions / insertions (Unreleased).** A segment-editing
        mod on an incoming PSF can't be layered onto a readpsf'd (immutable) segment, so `psfgen.do()`
        routes it to build-mode **full rebuild**: reconstruct the `Molecule` from the incoming PSF+coords
        *with* the mods (the build path's seqmod pipeline already applies mutations/deletions/insertions)
        and re-segment. A re-buildability guard hard-errors up front, naming any residue whose CHARMM
        `RESI` the release lacks (chosen over selectively preserving custom-residue segments); the
        incoming box (xsc) carries forward. Preserves sequence/coords and re-derives standard
        links/ssbonds; re-derives topology from RESIs (non-standard connectivity not carried).
        Validated on BPTI: `A:ASN,24,ALA` re-segments (disulfides re-derived), `A:56-57` deletes two
        residues, a bogus resname hard-errors. Replace/extend grafts now also work via this path
        (a small follow-on to route them there instead of the P2.5 error).
- [ ] **Special-position dedup is segment-granular.** Generating a biological assembly no longer
      replicates a segment that its transform maps onto itself — the fix for insulin 2ins, whose two
      axial zincs (deposited once at 1/3 occupancy on the 3-fold) were tripled into six, adding four
      spurious Zn(2+) and +8 charge. The invariance test is per *segment*, though, so a segment mixing
      on-axis and off-axis atoms is not invariant as a whole and is still replicated. In practice
      special-position atoms sit in their own segment (as in 2ins), so this has not bitten; closing it
      means testing per atom and filtering the authored PDB rather than skipping the whole segment.
- [ ] **Decide where `scripts/pestifer-snapshot` belongs.** A headless VMD renderer added for the docs
      figures: it picks representations from what a system contains (cartoon protein by chain or by
      secondary structure, licorice glycans, lipids with head-group beads, bulk solvent separated from
      ligands by residue copy count) and offers `-highlight`, `-face`, `-side`, `-ghost` and `-domains`
      for composing a specific figure. It lives in `scripts/` and so is repo-only — not installed by
      pip. Options: leave it as a dev tool, ship it as a `pestifer snapshot` subcommand, or fold the
      representation logic into a library module the subcommand and the docs build both call. Note the
      hard-won constraint recorded in its header: the Tcl must reach VMD on **stdin**, never `-e`, or
      representations are never built and every render silently produces an image containing nothing
      but the corner axes.
- [ ] **Drop the base.yaml dict-assertion test once ycleptic's `check-spec` ships a shape check.**
      `tests/unit/test_tasks/test_validate.py`'s schema walk, and the ad-hoc assertion that every
      `attributes:` entry is a dict, exist because `yclept check-spec` validated *vocabulary* and not
      *structure*: it passed, exit 0, on a `base.yaml` where an attribute indented four spaces too
      deep had been swallowed into its neighbour's `default:` list, so `charmmff.custom.prm` silently
      ceased to exist and `charmmff.custom.str` gained a dict where a filename belongs.
      ycleptic-repo fixed it in `cameronabrams/ycleptic#7` (branch `fix/speccheck-structure`), which
      detects a `default:` containing something shaped like an attribute spec. Verified here
      2026-10-01 against that branch in an isolated venv, gating on file content rather than the
      version string: pestifer's real `base.yaml` exits 0 (no false positive) and the slipped repro
      exits 1, naming the cause. **When it releases**, raise the `ycleptic>=` floor in
      `pyproject.toml` and delete the local workaround — but only after confirming the new
      `check-spec` is what actually ran, since a floor already satisfied by the installed version is
      precisely how a dependency test ends up testing PyPI instead (see CLAUDE.md, row 1).
      Report: `~/.local/state/fleet/drops/ycleptic-check-spec-shape-blindness-20261001.md`.
      **DONE 2026-10-01** (`75bec92b`), with one deviation worth keeping: the local assertion was
      not simply deleted, because nothing in pestifer *ran* check-spec — it was a line in
      `contributing.rst` asking contributors to remember a command. The manual step became
      `test_validate.py::TestCheckSpecCatchesAMisplacedAttribute` first.
- [x] **Raise the ycleptic floor again when the top-level element check ships (optional, small).**
      **DONE 2026-10-01**: ycleptic 2.4.4 released; floor raised and a second negative control
      added for the root list, proven to fail on 2.4.3 and pass on 2.4.4. Each control now fails
      below the floor it justifies.  Original note follows.
      2.4.3's structure check is applied to *nested* `attributes:` lists but not to the one at the
      root: `check_base_spec` iterates `base['attributes']` under `if isinstance(node, dict)` with
      no `else`, so a bare scalar at the root is silently ignored — the same defect class 2.4.3
      closed, one level up. Found here 2026-10-01 while testing the negative control that
      ycleptic's own docs recommend, which did not fire for exactly this reason. Fixed in
      cameronabrams/ycleptic#8 (`fix/speccheck-top-level-elements`), unreleased as of 2026-10-01.
      Verified against that branch: pestifer's real `base.yaml` still exits 0, so **no gate moves
      and this is not urgent** — bump only to close the narrow hole that a malformed root element
      in our own schema would currently slip through. pestifer's own control re-introduces a
      *nested* swallow and fires on 2.4.3 and on #8 alike, so it needs no change either way.
      Report: `~/.local/state/fleet/drops/ycleptic-ce30989-control-does-not-fire-20261001.md`.
- [ ] **The membrane-equilibrate convergence gate is not a convergence test — and the metric lives
      at TWO call sites.** Reported 2026-10-02 by pestifer-sweep with trajectory evidence:
      `~/.local/state/fleet/drops/pestifer-membrane-equilibrate-defects-20261002.md`, evidence trees
      under `~/devtests/pestifer/ex17-equil-analysis/`. **Not started — the replacement criterion is
      Cameron's design decision, not a cleanup.**

      Verified here from source, independent of their trajectories:

      - `_area_plateau_drift` (`tasks/membrane_equilibrate.py:231`) returns `(m2-m1)/m1` over the
        3rd vs 4th quarter of the stage-2 series. Those centroids are `T/4` apart, so for a steady
        descent the value scales as `slope x T` — **it grows with the window being tested**. A
        convergence test must get harder to pass as the run lengthens; this gets easier, then
        harder, then easier, which is why they measure 2–4 disjoint pass windows and why every
        gate-off run passes at the earliest evaluable point (`tail.size < 8`, i.e. 16 samples).
      - **The same arithmetic is implemented a second time**, independently, in
        `make_membrane_system._area_convergence` (`:859`) — `tail = x[len//2:]`, `h = tail.size//2`,
        `(m2-m1)/m1`. Its own docstring says the two "agree". They do: they are equally wrong.
        `_plateau_drift_for_report` falls back to it, so **fixing only `membrane_equilibrate` leaves
        the calibration-reliability check still using the broken metric.** This is the fourth
        instance of CLAUDE.md's "one idea at several call sites" pattern; sweep with
        `grep -rn 'tail.size // 2'` and treat the second site as part of the same fix.
      - The gate's observables are exactly `density` and `area` (`JointConvergence` at `:157`/`:163`)
        — no thickness, no chain order. Their Defect 4 stands.
      - `SAPL` is one `float` for the whole bilayer (`schema/base.yaml:1352`) and there is **no
        thickness key anywhere in the schema**, so their Defect 2 is a code fix by construction: a
        user whose bilayer builds 8 Å too thin has no way to say so.

      Resting on their trajectories, NOT verified here (I have no runs): that ex17's published
      47.06 APL is an artifact of where this metric first dips below tolerance from SAPL 50; the
      3–8 Å initial-thickness errors; patchB's stage-1 certification at rho=1.0101 while ~8 Å thin.
      Reference values for testing a fix are in their file. One correction to it: `compare.CHAIN_C`
      is in *their* analysis script, not in pestifer — pestifer ships no chain-atom table at all,
      which strengthens rather than weakens their argument for P–P thickness (phosphorus z needs no
      per-lipid table) over order parameter S.

      Their two load-bearing cautions, worth keeping attached: **ex16's quilt gate is correct**, so a
      replacement must not simply run everything longer; and **do not ship a seed fix alone** — it
      improves the numbers while leaving them seed-dependent, hiding the defect.
- [ ] **ex17's `lower_leaflet_phase: Ld` is a config defect — but fixing it recovers ~7%, not the
      deficit.** Measured 2026-10-02 with the build-time P-P logging from `46488ed2`.
      **Nothing changed in ex17; this is Cameron's call, and the sweep is mid-run.**

      ex17's config declares the lower leaflet Ld and justifies it: *"the cis double bond can't be
      ordered by the trans bias (its Lo ensemble == its fluid Ld one), so declare it Ld: same
      result"*. pestifer-sweep noticed SOPE/SOPS are stearoyl and POPE palmitoyl — **sn-1 is fully
      saturated in all three**; the cis kink is sn-2 only — so the trans bias should still extend
      sn-1 and the two ensembles should differ. Built patchB's composition both ways at ex17's
      SAPL 50:

          Ld   28.083 A      Lo   29.008 A      Lo - Ld = +0.924 A

      Identical across three seeds each (sd 0.000): `built_pp_thickness` is deterministic, since
      every anchor is pinned to one plane and the jitter is xy-only. So +0.924 A is an exact
      difference, not a sample. **The comment's "==" is false.**

      **And the code already said so.** `make_membrane_system._warn_phase_composition_mismatch`
      fires on this exact leaflet — 43% sterol, declared Ld — with *"expect the leaflet to order
      (drift toward Lo) during equilibration; consider phase: Lo"*. That prediction is precisely
      what the sweep measured (P-P 33.89 -> 42.30, still rising; corr(S, P-P) = +0.955). A correct
      warning was overridden by a comment whose premise is wrong.

      **But do not read this as the fix.** The grid builds patchB at 28.08 A against an
      equilibrated 42.3+; switching to Lo recovers 0.92 of ~14 A. Both ensembles are too curled,
      so the deficit is in how extended the MC conformers are, not in which ensemble is selected.
      Flipping ex17 to Lo would improve the number slightly, leave the defect, and spend the
      example's comparability to do it — the seed-fix-alone trap the sweep warned about, in a new
      costume.

      Open: (a) correct ex17's comment (false as written) — touches an example file the running
      sweep hashes, so not done unilaterally; (b) whether to flip the setting at all; (c) the real
      lever, conformer extension in the MC sampler. Probe kept at `~/devtests/pestifer/ld-vs-lo/`.
- [ ] **ROOT CAUSE: the MC conformer sampler confines footprint but not direction — 90% of drawn
      phospholipid conformers have a tail pointing out of the membrane.** Found 2026-10-02 chasing
      Defect 2. **Not fixed — changing the sampler changes every cached conformer and every
      membrane build, so it is Cameron's call.** Evidence:
      `~/.local/state/fleet/drops/pestifer-conformer-fold-20261002.md`, probes in
      `~/devtests/pestifer/ld-vs-lo/`.

      `athermal_mc.build_lipid_mc` confines tails to a **cylinder** about the membrane normal,
      which caps the in-plane footprint near the target APL. A cylinder is infinite along z, and
      nothing constrains a chain to the tail side of the headgroup, so a 180-degree fold back up
      is radially compliant, sterically allowed (self-avoidance only), and accepted. The sampler
      therefore hits its APL target honestly — which is why the conformer work was judged a
      success — while producing lipids with one tail in the water.

          phase Ld   POPE 10/10 folded   SOPS 8/10   SOPE 8/10
          phase Lo   POPE 10/10          SOPS 10/10  SOPE 10/10

      One POPE conformer, phosphate at z=8.08: sn-2 oleoyl runs C22:5.0 -> C218:17.1, exiting into
      the water 9 A above the phosphate, while sn-1 palmitoyl runs to -7.3 correctly.

      This is the thickness defect. Placement takes `raw_off = mean(anchor_z - mean(tail-tip z))`,
      which one folded tip collapses — POPE reads 3.30 A against 14.54 A measured to its deepest
      tip — and P-P follows directly, hence patchB's grid at 28.08 A against 42.3+ equilibrated.
      It also explains why declaring Lo bought only +0.92 A: the trans bias genuinely works (POPE's
      deepest tip goes 14.54 -> 21.65) but the other chain folds up more consistently, so the gain
      never reaches `raw_off`. The ordered ensemble is MORE folded than the fluid one.

      Note the regression was recorded as physics: `test_bilayer_spec_out`'s comment reads "(78.14
      with the athermal-MC conformer set, whose melted tails give a smaller z-extent than the old
      extended-rod conformers that gave 87.86.)" The z-extent collapse was read as melting; it is
      at least partly folding.

      **Refined 2026-10-02, second pass.** It is not sterol and not chain length. Fold rate:
      DMPC/DPPC/DSPC/POPC/PSM all 0/10; POPE 10/10, SOPE 8/10, SOPS 8/10. **POPC vs POPE is a
      controlled comparison** — identical chains, identical head/tail annotation, identical
      cylinder (`auto_cylinder_apl` sizes from chain count alone) — and one folds while the other
      does not. PE and PS fold; PC and SM do not. Why the headgroup matters is NOT established.

      Two layers permit it, and both must be fixed: the cylinder is infinite along z so nothing
      FORBIDS the fold, and the ordering objective `s = ½(3cos²θ − 1)` depends on cos²θ, so it is
      invariant under θ→180−θ and nothing PENALISES it — a folded chain scores exactly as well as
      an extended one. That predicts the observed Lo/Ld asymmetry (pushing the order target harder
      gives more folding, more reproducibly: Ld 26/30 sd 0.76, Lo 30/30 sd 0.39). It also closes a
      loop: `chain_order_parameter` is only valid for a head-up conformer, and the head-up
      canonicalization keys on the same tail-tip mean that a hairpin destroys.

      **Third pass: the fold is about HALF the deficit, not all of it.** Raw grids at ex17's
      SAPL 50 against gate-off equilibria: patchA (PSM/POPC/CHL1, no folders, Lo) 39.87 vs 46.48 =
      **+6.61**; patchA as Ld 36.48 = +10.00; patchB (PE/PS folders, Ld) 28.08 vs 41.12 = **+13.04**.
      So the headgroup/fold story explains the severity ORDERING, but a fold-free composition is
      still built 6.6 A thin. There is a second, fold-independent deficit — conformers are too
      curled even when no chain inverts, which is the "melted tails" z-extent collapse the
      `test_bilayer_spec_out` comment records. **Fixing the fold alone would take patchB from
      ~13 A to ~6.6 A short, not to correct.** Both need addressing.

      Corollary worth keeping: declaring Lo buys +3.39 A on patchA (non-folders) against +0.92 A
      on patchB (folders). The trans bias works as designed when chains cannot satisfy it by
      inverting — which is the cos²θ degeneracy showing up as a dose-response.

      **SEQUENCING — do not ship the P-P gate before the sampler fix.** Every membrane build
      currently starts 6.6-13 A thin, so a thickness gate would correctly refuse to certify almost
      all of them and convert every membrane build into a ceiling-length run (patchA took ~2.5M
      steps to climb 39.87 -> 46.48). ex16 is the converse existence proof: it builds near
      equilibrium, settles inside 180k steps, and **the EXISTING gate already gets the right answer
      there** — gated 59.37 ± 0.27 against gate-off 59.14.

      So the conformer fix may demote the gate change from "the fix" to "a safety net". That bears
      directly on the 2026-10-01 decision to forfeit 3.22.8 provenance in order to get a gate
      change: if the sampler fix lands first, the existing gate may be adequate for builds that
      start where they belong, and the forfeit may not have bought what it was spent on. Worth
      re-deciding before any provenance is actually given up.

      This is "never ship a seed fix alone", inverted: do not ship the detector without fixing the
      cause. (pestifer-sweep, `criterion/PP-REFERENCE.md` final section.)

      Fixes in order of confidence: (1) reject sampler moves putting a tail heavy atom above the
      head reference — the actual defect; (2) make `_mean_anchor_offset` robust to one folded chain
      (defence in depth, but alone it hides the sampler behind a better number); (3) reject folded
      conformers at load with a loud warning. **Unchecked:** ex16 is pure DMPC and its gate works;
      its conformers have not been examined and should be before generalizing.
- [ ] **Five lipids still build folded, because their IC-built STARTING structure is folded.**
      Found 2026-10-02 running down the residue after the axial penalty (`2e40b013`). Collection
      and evidence: `~/devtests/pestifer/conformers-regen/` (README has the full table).

      After regeneration the fold rate is 471/2440 (19.3%, 53 entries) -> 35/2440 (1.4%, 5
      entries). The five are `OSM`, `23SM__Lo`, `ASM__Lo`, `LSM__Lo`, `PMCL1`.

      They are **not a sampler failure**. Each starts with a chain 13-20 A above its own
      headgroup, built that way from internal coordinates, and the MC inherits it. `DMPC` is the
      control that settles it: it also starts folded, but only by 7 A, and the penalty pulls it
      fully out (-6.98 -> +18.43). Depth of the initial fold decides the outcome, and the penalty
      is visibly working on four of the five (23SM__Lo -20.19 -> -9.40) without finishing.

      **Raising `_AXIAL_BIAS` is the wrong lever and makes it worse**: 1.0 leaves 3 of 5 at 10/10
      with 2 partial; 5.0 and 20.0 put ALL five at 10/10 with worse margins. A large penalty is a
      wall, and a wall freezes a conformer in the basin it started in — the same mechanism that
      made the hard ceiling worse than useless on SOPE. So the fix belongs in the starting
      structure (unfold a chain built above the head plane, or try `--take-ic-from` /
      `--refic-idx`), not in the bias.

      Also corrected here: an earlier count of 16 affected entries, nine of them cardiolipins, was
      an **artifact of the audit script**, which looked for an atom named exactly `P` while
      cardiolipins name theirs `P1`/`P3`. Eight of nine cardiolipins are fine; only `PMCL1` really
      folds. Shipped code was never affected — `membrane_pp_thickness` selects phosphorus by
      element mass and excludes multi-phosphorus lipids explicitly.
- [ ] **Conformer duplication in the MC sampler (pre-existing, not from the axial work).** POPE
      yields 4 distinct `head-tail-length` values out of 10 under the old sampler and 5 of 10
      under the new one, so decorrelation is not delivering ten distinct conformers. Ensemble
      diversity is the whole reason the grid packer draws per-lipid rather than stamping one
      shape, so this is worth its own look. Noted while checking it was not a regression.
- [ ] **The on-demand autocache silently outranks an explicit `user_pdbcollections` override.**
      Found 2026-10-02 measuring the regenerated collection end-to-end.
      `PDBRepository.checkout` walks `self.registration_order[::-1]` — **reversed** — and
      provisioning registers base, then the user override, then the autocache dirs. So the
      autocache wins over a collection the user explicitly asked for. Both docstrings say the
      opposite: `pdbrepository.py:265` ("user-declared collections are searched first ... the base
      collection is searched last") and `:641` ("in the order they were registered").

      It is the silent-wrong-answer shape again: a user overriding a lipid that happens to be
      cached gets the cached one, and the build looks correct. It cost a measurement here —
      patchA read identically in both arms because its three Lo entries were cached, while
      patchB's were not and the override worked.

      **Not fixed:** reversing the precedence is a behavioural contract (an on-demand entry
      generated this run arguably should beat a stale user path), and it is Cameron's call which
      way it should go. Whatever is decided, the two docstrings and the code must be made to agree
      — right now no reading of them is correct.
- [ ] **patchA's residual thinness is NOT curl — it tracks sterol fraction, through the
      composition-weighted head plane.** Measured 2026-10-02 after the axial fix closed patchB
      (deficit 13.04 -> 2.11 A) while leaving patchA at 6.95 A.

      Curl is refuted, not assumed. patchA's chains sit at **88.9% / 88.5% of the all-trans
      ideal** z-rise, and its Lo ensembles reach their order target (PSM__Lo 0.297, POPC__Lo
      0.315, against a documented Lo band of 0.3-0.4). The conformers are extended and ordered.

      What pulls the plane down is the averaging. Placement sets
      `raw_off = composition-weighted mean over ALL species of (anchor_z - mean tail-tip z)` and
      pins every species' anchor to one plane. A sterol's anchor is its 3-OH, and it is far
      shorter than any phospholipid:

          PSM__Lo  19.35      POPC__Lo 22.93      CHL1__Lo 15.53      (A)

      patchA is 47% CHL1, so the sterol drags `raw_off` from 20.50 (phospholipids only) to 18.16
      -- predicting P-P 38.32 against 42.99, and reproducing the measured build of 39.53.

      Consistent across the three systems, ordered by sterol fraction:

          ex16    0% sterol   predicted 39.09, sweep measured build 39.91, relaxes DOWN (too THICK)
          patchB 43% sterol   1.06 A/leaflet short after the axial fix
          patchA 47% sterol   3.48 A/leaflet short

      **This is the cholesterol correlation pestifer-sweep originally reported**, which I retired
      as a confound when the FOLD turned out to be headgroup-driven. Both were right, for two
      different mechanisms: PE/PS folding drove patchB's catastrophic case, sterol fraction drives
      the fold-free residue and ex16's correctness.

      **Caveat before anyone acts.** The equilibrated reference already contains cholesterol's real
      effect on thickness, so "the build is thin" is measured against that. Whether co-planar
      anchoring of a sterol 3-OH with a phosphate is actually WRONG needs the equilibrated
      structures: in a real bilayer the 3-OH sits ~3-6 A below the phosphate plane. pestifer-sweep
      holds those trajectories; I do not. Ask before changing placement semantics.
- [ ] _(add items here)_
