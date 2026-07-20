# Changelog

Pestifer follows [Semantic Versioning](https://semver.org/) and documents changes below.

## [Unreleased]

- fix: **glycan coordinates no longer get silently scrambled during a build.** psfgen's `writepdb` emits CHARMM residue names verbatim, and six-character carbohydrate names (`BGLCNA`, `ANE5AC`, ...) overflow the four-column PDB `resName` field, shifting the x/y/z columns a few characters to the right. Every fixed-column reader downstream — VMD's own PDB reader *and* pestifer's numpy coordinate readers (`coorddf_from_pdb`) — then read the wrong substring, so the first pipeline step to load and rewrite that PDB (a crotation, a declash, or the graft-time ring-piercing resolver) destroyed the glycan coordinates (bonds stretched to hundreds of ångström). The psfgen scripter now re-anchors every written PDB's coordinate columns into standard positions (new `util/coord.standardize_pdb_columns`, applied in `PsfgenScripter.runscript`); the identity/occupancy/beta/segname columns are untouched and already-standard records are left byte-for-byte identical (idempotent). Verified on the 6CM3 BG505-SOSIP+sCD4 build: 1344 broken glycan bonds → 0, max glycan bond length 240.9 Å → 1.5 Å.

## [3.8.0] - 2026-07-17

- change: **all mmCIF parsing is offloaded to pidibble; the `mmcif` and `mmcif-pdbx` dependencies are dropped.** pestifer had carried two structure parsers — pidibble for PDB, and `mmcif`/`mmcif-pdbx` (via `util/cifutil.py`'s `CIFdict`/`CIFload`) for mmCIF — with each object's `from_pdb`/`from_cif` paths diverging. pidibble ≥1.7.1 normalizes mmCIF into the same PDB-record namespace, so the raw-CIF layer is deleted and pestifer now has **one** structure parser. **Requires `pidibble>=1.7.1`** (for `header_label` on REMARK.350 records). mmCIF label-primary chain identity is preserved exactly, including the `cif_residue_map` artifact. The one behavior change: mmCIF symmetry operators normalize `1_555` → `1555`, which *aligns* CIF with the PDB path (which already emitted `1555`); `sym1`/`sym2` are only ever written to SSBOND/LINK `pdb_line` output. Validated by a golden AsymmetricUnit oracle on 4zmj/6pti, the full unit suite, and an all-examples runthrough. See `docs/design/mmcif-offload.md`.
- fix: **examples 16 and 17 (membrane) build again.** `make_membrane_system`'s orientation step tested `orient.get('source') is not None`, but the schema materializes an empty `orient: {source: [], target: []}` default even when only the `z_head_group`/`z_tail_group` shorthand is given — so the empty `orient` branch was taken and the `ALIGN` operation failed schema validation before the protein was ever oriented. The check is now for a non-empty source/target, letting the shorthand take effect.
- fix: **example 25 (multi-step continuation build) runs standalone.** `build-example` now runs an example's `inputs/aux/*.yaml` helper configs first, in sorted (dependency) order, so their outputs (`base.psf`, `fab_p*.psf`, ...) exist before the main script's `continuation`/`merge` consumes them. Previously `build-example 25` failed immediately with `No such file or directory: 'base.psf'` because the helpers had to be run by hand. No effect on examples without an `aux/` folder.
- enhancement: **pestifer warns when an available GPU is going unused.** GPU-resident NAMD is elected only when `paths.namd3gpu` differs from `paths.namd3`; both default to `namd3`, so on a host whose GPU-resident build is a *separate* binary (e.g. `namd3gpu` alongside a multicore `namd3`) every MD task silently ran CPU-only. pestifer now warns — naming the candidate binary and the setting to change — when a GPU is detected, a GPU-resident build is on PATH, and `paths.namd3gpu` still points at the CPU binary. Silent when there is no GPU, no separate GPU build, or when both paths resolve to the same binary.
- docs: new **Example 26 — Subtilisin Carlsberg in acetonitrile**, completing the organic-cosolvent series with Examples 23 (DMSO) and 24 (acetone). Like acetone, the acetonitrile box is generated on demand; it also exercises the degenerate-torsion fill (acetonitrile's linear H–C–C≡N dihedral, supplied as `k = 0`) so the on-demand box equilibrates.
- enhancement: **on-demand solvent boxes provision parameters robustly, and acetonitrile (`ACN`) now builds.** When `make_solvent_box` equilibrates a freshly-packed box, it now (A) also loads a residue's companion parameter file if the residue is defined in a bare `.rtf` (a `.str` is self-contained), and (B) turns a missing-parameter NAMD abort into an actionable error naming the exact residue and atom-type term — instead of a cryptic `UNABLE TO FIND ... PARAMETERS` fatal mid-equilibration. Genuinely-degenerate torsions that CGenFF omits but psfgen autogenerates — e.g. acetonitrile's `HGA3 CG331 CG1N1 NG1T1`, degenerate because the C–C≡N skeleton is linear — are covered by a curated, reviewed `k=0` entry in the new `custom/toppar_pestifer_dihedral_fills.prm`, so `ACN` boxes out of the box. This is *not* a blanket zero-fill: any other unrecognized missing term is still surfaced as an error, never silently masked. (Discovered building Example 24; see `docs/design/on-demand-solvent-parameters.md`.)
- change: **the bundled CHARMM force field is now a single release, `February2026` (again the default), and `July2024` is removed.** 3.7.2 had temporarily shipped only `July2024`; this consolidates on the newer `February2026` (CGenFF v5 superset) as the sole bundled release. **Action required if you pinned it:** a config with `charmmff.release: July2024` now hard-errors — remove the pin to use `February2026`. Unpinned builds move from the July-2024 back to the February-2026 force field. (If you have a stale local parsed-FF cache, run `pestifer cache rebuild`.)
- build: **the bundled toppar tarball is trimmed to what pestifer actually loads.** New `scripts/trim_toppar.sh` strips an upstream `toppar_c36_<release>.tgz` down to the top-level topology/parameter/stream files plus the `stream/` subtree, dropping the subdirectories pestifer never reads (`metals/`, `drude/`, `non_charmm/`, `silicates/`, implicit-solvent sets, etc.). Verified to produce byte-identical build topology *and* identical residue/patch/collection counts — the dropped files were never indexed. Trims feb26's toppar from 14 MB to 5.3 MB (jul24 was 11→3.1 MB); with `July2024` also removed, the wheel/sdist drop from ~15 MB to ~12 MB.

## [3.7.2] - 2026-07-14

- change: **the default CHARMM force-field release is now `July2024`, and `February2026` is no longer
  shipped** in the package (kept in the source tree, excluded from the build). This roughly halves the
  bundled force-field payload. **Action required if you pinned it:** a config with
  `charmmff.release: February2026` will now hard-error ("version directory not found") — remove the pin
  to use `July2024`, or re-add the `feb26` files. Builds that never set `charmmff.release` silently
  switch from the Feb-2026 force field to July-2024; pin `charmmff.release: July2024` to lock it.
  (If you have a stale local parsed-FF cache, run `pestifer cache rebuild`.)
- build: **`tests/` and `docs/` are excluded from the built package.** The source distribution had been
  bundling the entire test tree — ~175 MB of large PSF/PDB fixtures — because the excludes only covered
  a subset; none of it is needed to install pestifer (the wheel already omitted it). Together with
  dropping `feb26`, the sdist shrinks from ~65 MB to ~15 MB (and the wheel likewise), keeping releases
  well under PyPI limits and slowing growth of the project's total-size quota.

## [3.7.1] - 2026-07-14

- new feature / change: the **`ANGLEIJK` rotation is re-homed as the `AXISANGLE` movetype of the `transrot` mod**. `ANGLEIJK` rotates a whole segment as a rigid body (about the axis normal to a three-atom angle), so it belongs with the rigid-body `transrot` operations rather than the internal-coordinate `irotations`. The new `transrot` `AXISANGLE` movetype takes `axis_atoms` (three atomselections `[selI, selJ, selK]`) + `angle` (+ optional `sel`), pivoting at `selJ` about `(rI-rJ)×(rJ-rK)`. `ANGLEIJK` remains accepted as an `irotations` entry for backward compatibility but is **deprecated** (it logs a warning pointing to `AXISANGLE`). The rotation geometry is now a single shared engine used by both paths.
- change: the **`crotations` mod is renamed `irotations`** ("i" for *internal*-coordinate rotations), in both the `manipulate` and `psfgen` mod sets. **`crotations` is retained as a backward-compatible alias** (mirroring the `transrot`/`rottrans` pattern), so existing configuration files keep working unchanged — either spelling validates and dispatches to the same handler. The shortcodes (PHI/PSI/OMEGA/CHI1/CHI2/ANGLEIJK/ALPHA/GLYCAN_PENDANT) are unchanged. (Also fixes a stale `manipulate`-task key that never matched the mod's header.)
- change: pestifer's **built-in `custom/` CHARMM force-field files now live in a single shared `charmmff/custom/` directory** (a sibling of the per-release version directories) instead of being duplicated under each release. The files are release-independent, so this de-duplicates them; `patches/` remain per-release. `modify-package charmmff add-residue` now installs into the shared directory. (The `feb26` and `jul24` copies were identical apart from a few descriptive comments, which were retained.)
- change: the **`domainswap` task is retired** and no longer available in a build config. A config using `domainswap:` now fails validation (it is removed from the task schema and the task registry). The code is preserved in the source tree (`pestifer/tasks/domainswap.py`, kept importable but unregistered) and can be reinstated by re-adding its schema block and registry entry; its documentation pages are retained as orphaned reference.
- bugfix: **branched glycans had their branch-point bonds wired to the wrong residues.** A glycan segment stores its residues in BFS (root-outward) order and renumbers them in that order, so the `LINK`-derived psfgen patches reference BFS resids; but the generic-segment coordinate writer emitted `$sel set resid [list ...]` with the resid list in BFS-residue order while VMD's `atomselect` returns atoms in ascending-serial order and applies the list in *that* order. Past the first branch point the two orderings diverge, so the written coordinates were numbered inconsistently with the topology and the branch bonds spanned >10 Å (NAMD then failed on the impossible geometry). Serials are now paired with their resids and sorted by serial so the list matches VMD's atom order — a no-op for serial-ordered protein segments, a correction for branched glycans. Verified on 2wah/4b7i/4byh (added as an integration regression test) and 6cm3, whose ASN-234 high-mannose glycans now match their `LINK` records exactly. (Surfaced misbuilding the ASN-234 glycans of PDB 6cm3.)
- bugfix: **chloride ions were left at the origin.** The built-in `CL → CLA` residue alias builds the CHARMM `CLA` residue (whose atom is named `CLA`), but RCSB PDBs name the chloride atom `CL`; with no matching atom alias, psfgen could not map the coordinate and dumped the ion at `(0,0,0)`. Because `CL` was already residue-aliased, the monatomic-ion pre-build guard also skipped it. Added the missing `CLA CL CLA` atom alias so chlorides build in place. (Zinc needs none — the `ZN2` residue's atom is already named `ZN`.)
- bugfix: **`new-system` leaked structure-specific directives into every generated config.** The subcommand scaffolded from `self.examples[0]`, but the examples list is populated in filesystem-iteration order, so `[0]` was nondeterministic and happened to resolve to the subtilisin example — injecting its `CA CAL` / `CAL CA CAL` calcium aliases (and, with `--full`, its DMSO solvent and validate tests) into unrelated new systems. `new-system` now scaffolds from a dedicated, structure-neutral template (`resources/templates/new_system.yaml`).
- bugfix: **the `ligate` task duplicated the peptide bond of a built N-terminal loop**, producing a topology NAMD rejected with a degenerate `C NH1 C` angle (no parameters). When `build_zero_occupancy_N_termini` builds a missing N-terminal loop, psfgen bonds it *contiguously* to its first resolved neighbor (no terminus break, since the sacrificial-residue break is only made for interior loops), but the gap-healing pass (`write_gaps`/`write_connect_patches`) still emitted a heal `LINK` for it — re-creating the already-existing peptide bond. Healing is now restricted to interior loops (`0 < i < len-1`), matching where psfgen actually creates a break; interior loop healing is unchanged. (Surfaced building PDB 6cm3 with its gp41 fusion-peptide N-termini enabled.)

## [3.7.0] - 2026-07-10

- bugfix: `merge` produced a system in which **merged copies of the same structure collapsed onto one chain ID** in any chain-based view, so three completed 17b Fabs (each entering as chain H/L) displayed as a single Fab even though all their atoms and segments were present and correctly placed. The cause is that a PSF carries no chain column: every time the pipeline regenerates a PDB from the PSF (the `solvate` step, and every `coor`→`pdb` conversion), VMD re-derives each atom's chain ID from its segid's leading characters — so segids `H`, `H0`, `H1` all become chain `H`. `merge`'s segid-enumeration (`H` → `H0` → `H1`) therefore did not actually separate the copies downstream. `merge` now renames each colliding segment to a **fresh single-character segid** whose leading character is unused, so its derived chain ID is unique *and* stable through every regeneration (verified end-to-end: three Fabs survive solvation as three distinct chains). (Surfaced completing the three partially-resolved Fabs in Example 25.)
- enhancement: **multiple builds in one directory no longer clobber each other's diagnostics log.** `pestifer build`/`run` now derive the default log name from the input config — `build foo.yaml` writes `foo-diagnostics.log` instead of a shared `pestifer_diagnostics.log` — so successive builds keep distinct logs automatically (a global `--log-file` still overrides). Other subcommands keep the `pestifer_diagnostics.log` default; the two subcommands that already set their own (`density-profile`, `mdplot`) are unchanged. The misleading `--log-file` help ("default: subcommand-specific") is corrected.
- enhancement: for the same reason, the `terminate` task's **build-intermediate artifacts tarball now defaults to `{basename}-artifacts.tar.gz`** (from the final-system `basename`) instead of a fixed `artifacts.tar.gz`, so multiple builds in one directory produce uniquely-named tarballs. An explicit `artifacts:` still wins, and a build with no basename falls back to `artifacts.tar.gz`.

- bugfix: an example with **multiple YAML input scripts** (a main script plus helper scripts, as in Example 25) was silently **dropped from the example registry** — the scanner keyed off the sole YAML in an example's `inputs/` directory and skipped any directory with more than one, so `fetch-example`/`run-example` reported it "not found". Auxiliary inputs now live in an `inputs/aux/` subfolder, leaving exactly one main script in `inputs/`; the example is discovered, and `fetch-example` retrieves the main script *and* all `aux/` helpers into the working directory.

- docs: new **Example 25 — completing partially-resolved ligands**. Builds the full sCD4/17b-liganded HIV-1 Env trimer from 5vn3, whose three 17b Fabs are only half-resolved (variable domains only). A tutorial-style, multi-script workflow completes each Fab from a fully-resolved second structure (1gc1) — `align` the complete Fab onto the partial Fv, `transfer_coords` to overwrite the Fv with the exact 5vn3 interface coordinates, then `merge` the three completed Fabs into the Env+sCD4 base. Demonstrates cross-structure alignment/coordinate-transfer with mismatched residue numbering, and automatic segment-name collision handling in `merge`. (Companion to Example 12, which builds the Env alone.)

- bugfix: the `merge` task produced a PSF that **NAMD rejected with `DIDN'T FIND "NATOM"`** whenever it injected or stripped header REMARKS — i.e. any merge of disulfide-bearing systems whose segments were renamed to resolve a collision (the renamed patch remarks are re-injected, and dropped topology remarks are removed). Both edits changed the number of title lines without updating the `!NTITLE` count, so NAMD read the wrong number of title records and never reached `!NATOM`. `_inject_patch_remarks` and `_strip_topology_remarks` now recompute `!NTITLE` after editing. (Surfaced merging three completed antibody Fabs, each carrying disulfide patches, into an Env trimer.)
- bugfix: a `manipulate` **`align` or `transfer_coords` selection-count mismatch is now a hard error** instead of a silent failure. The Tcl congruency check `exit 1`s on a mismatch, but VMD returns process code 0, so the task previously reported success while writing no output pdb — the failure only surfaced at the next task as a cryptic "cannot open file". The congruency errors now carry the `PESTIFER-ERROR` marker that the manipulate task already scans for, so a mismatch aborts the build immediately with the offending selections and their atom counts named.

## [3.6.0] - 2026-07-09

- enhancement: `show-resources resname` now displays the residue's **long name** — the free-text comment after `!` on the CHARMM `RESI`/`PRES` title line, typically a molecular formula and/or descriptive name (e.g. `BORO` → "B1O2C1H5, methyl boronic acid, neutral"). It appears as a quoted line under the topology entry in the detailed view and at the end of each `--contains` search-result line, making the search far more informative. (Sourced from the parsed force field, so it takes effect on the next release's cache rebuild or via `pestifer cache rebuild`.)
- change: the bundled **February 2026 CHARMM force field ships two CGenFF versions** (v5 and the older v4.6) as parallel `top_all36_cgenff*.rtf` / `par_all36_cgenff*.prm` files; pestifer now loads **only v5** and skips the `_v4.6` files. v5 is a strict superset — every v4.6 residue and atom type is present in v5 (the handful of shared residues that differ do so only in comments/whitespace) — so loading v4.6 merely duplicated ~937 residue definitions. (Takes effect automatically on upgrade, since the parsed-force-field cache is keyed by pestifer's major.minor version; or run `pestifer cache rebuild`.)
- new feature: `psfgen` now **detects a common monatomic-ion residue-name collision before it builds** and stops with the exact fix. A bare (one-atom) input residue named `CA`, `NA`, `K`, `CL`, `ZN`, `CS`, `RB`, `LI`, or `CD` is unambiguously a metal/halide ion, but psfgen would otherwise build the same-named CHARMM/CGenFF residue (for `CA`, a 68-atom organic molecule) and misplace its atoms. When such an ion is present and not already aliased, the task aborts *before* the (wasted) build with a message naming the residues and the residue+atom aliases to add (e.g. calcium → `residue: ["CA CAL"]`, `atom: ["CAL CA CAL"]`). It never silently reinterprets your input, and it only fires for genuinely bare ions (a real multi-atom residue of the same name is left alone). The origin-coordinate guard below remains the general backstop for any unrecognized cause.
- bugfix: a `psfgen` build that leaves atoms at the **origin (0,0,0)** is now a **hard error** instead of a silently-accepted broken structure. psfgen writes any atom whose coordinates it can neither read nor guess at exactly `0.000 0.000 0.000`; this classically happens when a PDB residue *name* collides with a different, multi-atom CHARMM/CGenFF residue of the same name — most often a monatomic metal ion (e.g. calcium named `CA`, which matches a CGenFF organic residue), so psfgen builds the organic residue's extra atoms with no coordinates for them. Previously the pipeline treated this as a successful build and the following minimization blew up. The task now scans the psfgen output PDB, and on finding origin atoms aborts with a message that names the offending residues and shows the alias fix (for calcium, `residue: ["CA CAL"]` / `atom: ["CAL CA CAL"]`). Note psfgen's own "poorly guessed coordinates" warning does *not* flag these atoms, so the coordinate scan is what catches it. (Surfaced building Example 23 — subtilisin (1scd) has two calcium ions.)
- docs: new **Example 24 — Subtilisin Carlsberg in acetone**, a companion to Example 23 (subtilisin in DMSO) that showcases **on-demand solvent-box generation**: acetone (``ACO``) is not a shipped solvent, so the ``solvate`` task builds and caches its box on the fly. Built end-to-end (auto-generated box: edge 30.0 Å, ρ 0.77; system neutralized by replacing 3 acetones with Cl⁻).
- bugfix: ``modify-package example add`` given a YAML **path** (rather than a bare ``<name>.yaml`` in the CWD) wrote a corrupt toctree entry and resource key — the full path leaked in as the example ``shortname``. The shortname is now always the bare basename, and the source path is threaded to the file-copy step, so a path outside the CWD works.
- bugfix: **on-demand solvent-box generation** (v3.5.0) never actually ran in a real `solvate` build. The `solvate` task provisions only the PDB repository (not the force field's residue objects), so the miss handler's `if solvent not in CC` membership check raised `CHARMMFFContent must be provisioned before accessing residues or patches` before generation could start (the feature's unit tests mocked the generator and its pre-flight used a fully-provisioned content object, so both missed it). The force-field-membership check now happens inside `ensure_solvent_box`, which builds in a fresh, fully-provisioned `ResourceManager`; a non-force-field solvent surfaces as a clear error there. Verified end-to-end by building a subtilisin system in acetone (a non-shipped CGenFF solvent), which now auto-builds and caches the box. (The lipid-conformer path was unaffected — the membrane packer provisions its content before looking up species.)
- new feature: the `psfgen` task now **resolves glycan ring-piercings at graft time**. A grafted or model-built glycan bond can thread through the hole of a protein or glycan ring with no atomic overlap; such a topological piercing has a zero clash count, so the existing clash-count glycan declash is blind to it, yet it survives minimization and RATTLE-fails at the first dynamics step. After building glycans, pestifer now scans for these piercings (winding-number test, reusing `RingChecker`) and clears each by rotating the offending glycan sub-branch about an upstream rotatable bond — the glycan is moved, the protein is not — before the structure is written. This runs whenever the system has glycans, independent of the declash `maxcycles`, and is best-effort/non-fatal (a piercing no rotation clears is warned, leaving a downstream `ring_check` to attempt further resolution, e.g. aromatic side-chain rotation). It makes a separate glycan `ring_check` unnecessary in the common case. Opt out with `source.sequence.glycans.declash.check_piercings: false`. Validated on a real fully-glycosylated 4zmj model: two glycan piercings (a proline and an aromatic His ring, both speared by glycan bonds) are cleared to zero. The rotation engine is now shared (`pestifer/psfutil/ring_resolve.py`) between this path and the `ring_check` task

## [3.5.0] - 2026-07-09

- new feature: **on-demand generation cache** for missing PDB-repository coordinates. When a `solvate` task requests a non-water solvent that is defined in the CHARMM force field but has **no shipped box** (previously a hard error telling you to run `make-pdb-collection`), pestifer now **builds the box on the fly** with the existing `make_solvent_box` pipeline and caches it per-user, keyed by CHARMMFF release, under `~/.pestifer/pdbrepository/<release>/solvent/<RESI>/`. The first build is a one-time cost (pack + minimize + NPT equilibration, loudly logged); the cache is auto-registered as a user PDB collection at provisioning, so every subsequent build reuses it instantly. Generation is concurrency-safe (an `fcntl` advisory lock serializes builders — auto-released if one dies — and the finished entry is published with an atomic rename) and runs in an isolated `ResourceManager` so the live build's force-field state is untouched. Cached entries are marked `quality: auto` in their `info.yaml` to distinguish them from hand-curated shipped boxes. Opt out with the new `charmmff.generate_missing_coordinates: false` (default `true`), which restores the former explicit error for reproducibility-strict or offline runs. This complements — does not replace — the contribute flow (`make-pdb-collection` + `modify-package pdb-repo add-entry`), which is still how a curated box gets into the shipped package.
- new feature: the same on-demand generation cache now covers the **grid membrane packer**. A `make_membrane_system` composition that names a lipid (or other membrane component) defined in the CHARMM force field but absent from the shipped PDB repository — previously the hard error `Cannot find <RESI> in PDB repository` in `bilayer.py` — now **generates that residue's single-molecule conformers on the fly** with the existing `do_resi` sampler (`kind: molecule`) and caches them under `~/.pestifer/pdbrepository/<release>/lipid/<RESI>/`, reusing them on subsequent builds. The *kind* of artifact is driven by the consumer: the packer gets conformers, the `solvate` task gets a box. Same concurrency-safe atomic publish, isolated-builder, loud-one-time-cost, and `charmmff.generate_missing_coordinates` opt-out as the solvent-box path
- bugfix: `test_charmmffcontent_full_provisioning` asserted stale PDB-collection counts (lipid 152, solvent 12) that predated the shipped-collection growth (the 67-entry lipid batch-install and the MEOH/ETOH/DMSO solvent boxes); updated to the current 219/15

## [3.4.0] - 2026-07-09

- new feature: `make_membrane_system`'s embed step gains a general **`embed.orient`** spec that orients the protein via the same **`ALIGN`** operation as the `manipulate` `transrot` mod -- rotate a `source` vector onto a `target` vector (each a literal 3-vector or an atomselection pair) by the minimal roll-free rotation about the protein's center of mass. This lets a membrane build target something other than `+z` (e.g. a tilted normal) or use a literal source vector, all in the one task. The existing `z_head_group`/`z_tail_group` shorthand is preserved and now maps onto this path (`source: [z_tail, z_head], target: [0,0,1]`), so existing membrane configs are unchanged. Internally this replaces the bespoke `Orient::orient`-based `bilayer_orient` orientation with the shared, tested `transrot` `ALIGN` code path -- verified coordinate-identical (RMSD < 0.01 Å) on the van3 fixture, and additionally robust to the degenerate antiparallel case the old path did not handle

## [3.3.0] - 2026-07-09

- new feature: a ``manipulate`` ``transrot`` gains an **``ALIGN`` movetype** that rotates a fragment so a *source* vector is carried onto a *target* vector -- the vector-based orientation win for membrane embedding (orient a protein's spanning axis onto the membrane normal ``[0,0,1]`` without a reference structure). Each of ``source``/``target`` is either a literal 3-vector ``[x,y,z]`` or a pair of atomselections ``[selA, selB]`` whose mass-weighted centers define the vector. The applied rotation is the **minimal (roll-free) rotation** -- about the axis ``source x target`` by ``atan2(|source x target|, source . target)`` -- performed about the fragment's own center of mass, so no spurious twist about the target axis is introduced (degenerate already-aligned and antiparallel cases are handled). Reuses the ``sel`` selector and disconnection guard
- new feature: a ``manipulate`` ``transrot`` can now act on **a fragment of the system** instead of always the whole structure. Give the directive in dictionary form with a ``sel`` field holding a VMD atomselection (default ``all``) to translate/rotate just that fragment (rotations still about the fragment's own center of mass). Because a rigid-body transform only makes sense for a fully disconnected fragment, a **disconnection guard** is emitted for any selection other than ``all``: if any selected atom is bonded across the selection boundary the build **hard-errors** instead of silently stretching the crossing bond. (VMD's Tcl ``exit 1`` always returns process code 0, so the guard signals failure via a ``PESTIFER-ERROR`` marker that the ``manipulate`` task detects in the VMD log and turns into a non-zero task result.) The header ``rottrans`` is now accepted as a synonym for ``transrot``
- bugfix: a ``manipulate`` ``transrot`` ``ROT`` now rotates about the **center of mass** of the moving selection, as intended. The generated command used VMD's ``trans origin {COM}``, which *moves* the COM to the global origin and rotates there (so the structure was rotated about, and translated onto, the origin); it now uses ``trans center {COM}``, which rotates in place about the COM

## [3.2.1] - 2026-07-08

- bugfix: the ``manipulate`` task's ``transrot`` mod (rotations/translations of the whole structure) did nothing -- ``coormods`` dispatched on the object key ``rottrans``, but objects are keyed by their YAML header ``transrot``, so ``write_rottrans`` was never called and the generated script contained no rotation/translation operations. Fixed the key; also made ``write_rottrans`` and ``write_crot`` accept the ``molid`` argument that ``coormods`` passes (they previously would have raised ``TypeError`` once reached). Adds unit and integration regression tests
- bugfix: the ``task-table`` docs directive labeled the ``solvate`` task "water box" regardless of solvent; it now reports the actual solvent (e.g. "DMSO box") for a non-water build

## [3.2.0] - 2026-07-08

- bugfix: tiled non-water solvent molecules are **no longer malformed**. The self-contained topology pestifer generates for VMD ``solvate``'s ``-stop`` omitted the ``AUTOGENERATE ANGLES DIHEDRALS`` directive; since a CGenFF RESI lists only atoms and bonds and solvate's replica ``segment {residue <solvent>}`` sets no ``auto`` option, every tiled solvent molecule was built with bonds but **no angle or dihedral terms** and distorted freely under MD (a DMSO system equilibrated to ~1.27 g/cc instead of ~1.1). The generated topology now emits the directive, so replicas get their full angle/dihedral set (verified: a tiled DMSO has 15 angles / 12 dihedrals). The pre-equilibrated boxes themselves were unaffected
- bugfix: solvating with a non-water solvent **now fills the whole simulation cell** instead of placing a single box. Two problems: (1) VMD's ``solvate`` plugin hardwires the segid ``QQQ`` for the box it tiles, so a box built with any other segid was placed once and never replicated -- the built-in ``MEOH``/``ETOH``/``DMSO`` boxes (and ``make-pdb-collection solvent`` output) now use segid ``QQQ``; and (2) solvate builds its replica residues in an isolated psfgen context that needs the solvent's topology (RESI + atom-type masses), which a stock model-compound stream may omit -- pestifer now generates a small self-contained topology for the solvent and passes it via ``-stop``. Verified: a DMSO-solvated 1scd now yields ~1,800 DMSO in a ~320,000 Å³ cell (was ~210), and MEOH/ETOH fill correctly too
- new feature: pestifer can now **add ions to a non-water solvent** by replacing solvent molecules with ions -- a solvent-agnostic analogue of VMD's ``autoionize`` (which only works for water). The ``solvate`` task neutralizes the system with counter-ions by default (new ``neutralize`` spec, default ``true``); ``salt_con`` additionally adds cation/anion pairs sized from the box volume (default ions ``SOD``/``CLA``; multivalent ion valences are honored). For water this is autoionize as before; for a non-water solvent the ions are placed at the positions of solvent molecules chosen far from the solute and from one another (default 5 Å), which are deleted to make room, then an ``ION`` segment is built with psfgen (new ``PestiferIonize`` Tcl package). Set ``neutralize: false`` (and no salt) to keep the system's net charge instead. Validated: a DMSO-solvated +3 protein neutralizes to net 0 by replacing 3 DMSO with 3 ``CLA``, and runs through minimize + NVT
- bugfix: solvating with a **non-water solvent no longer fails at ionization**. VMD ``autoionize`` places ions by carving space out of *water*, so it cannot ionize a dense non-aqueous box (it errored "Failed to add ions after 10 tries" on a DMSO-solvated protein). Non-water solvents are now ionized by solvent replacement instead (see above)
- bugfix: the ``solvate`` task no longer **silently succeeds when it produced no output**. VMD exits 0 even when ``solvate``/``autoionize`` catch their own errors, so a failed ionization previously left the task registering a psf/pdb that were never written (the build "succeeded" with only an ``.xsc``). The task now verifies the expected psf/pdb exist and aborts with a clear message (pointing at the non-water-ionization cause when applicable)
- enhancement: the end-of-build **System Report** (``terminate`` task) now reports the **total system charge** (summed from the PSF per-atom charges) alongside the atom/bond/angle counts -- a quick check that a built system is net-neutral (or carries the expected net charge) before simulation
- new feature: ``modify-package`` now keeps a **modification ledger** -- an append-only, git-tracked record at ``pestifer/resources/modifications.jsonl`` (one JSON object per line). Every mutating command appends an entry (category, verb, one-line summary, files touched, committer, branch, unique id) that is committed alongside the change, so it travels in the contribution's branch/PR. It complements the ``show-resources`` provenance tags: provenance says a file is ``custom``, the ledger says who added it, when, and via which operation. New ``modify-package ledger`` category: ``ledger show [--limit N] [--category C]`` lists entries, and ``ledger revert <id>`` reverses a recorded modification -- it locates the commit that made the change, ``git``-reverts it into the working tree on a fresh branch, and curates the ledger (the original entry is marked reverted and a new ``revert`` entry appended) so the audit trail survives. Revert works uniformly for any verb and surfaces conflicts via git; only committed modifications (the default flow) are revertable this way
- change: ``modify-package`` now **opens a contribution branch and commits by default**, on every category including ``example`` (previously the git workflow was opt-in via ``--branch`` and only on ``pdb-repo``/``charmmff``). Each invocation verifies a clean working tree, makes the change on a fresh branch auto-named ``modpkg/<category>-<verb>-<detail>``, and commits exactly the files it touched. ``--branch NAME`` names the branch explicitly; the new ``--no-branch`` skips branching/committing and just applies the change to the working tree (the old default). Example operations, which don't enumerate their touched files, are committed via a ``git status`` diff (safe because the clean-tree precondition attributes every change to the operation, and the resource cache lives outside the repo). Adds ``gitutil.changed_paths`` and a branch-flow test suite
- bugfix: ``modify-package example add`` now works, and no longer corrupts the examples documentation. Three latent bugs: (1) the ``example add`` verb dispatched to a nonexistent ``ResourceManager.add_example`` (the method is ``append_example``), so it always raised ``AttributeError``; (2) ``gitutil.changed_paths`` parsed git-status output that had been ``.strip()``-ed, dropping the first character of the first path when it was a modified tracked file (``docs/...`` -> ``ocs/...``) and breaking the commit; (3) appending an example to the toctree hardcoded a two-space indent (mismatching the toctree's three-space options, so docutils read the options as entries) and deleted the blank line before the following section. All three are fixed and regression-tested
- new feature: the built-in ``solvent`` PDB collection now ships **pre-equilibrated boxes for three non-water solvents** -- methanol (``MEOH``), ethanol (``ETOH``), and DMSO (``DMSO``) -- so ``solvate: {solvent: MEOH}`` (and ETOH/DMSO) works out of the box with no ``make-pdb-collection solvent`` step. Each was built at nmol=216 with a full production NPT equilibration: MEOH ~24.6 Å / 0.77 g/cc, ETOH ~27.5 Å / 0.80 g/cc, DMSO ~29.5 Å / 1.09 g/cc. (Added to the newest bundled force field, ``feb26``.)
- bugfix: ``make-pdb-collection solvent``'s Open Babel coordinate fallback no longer adds spurious hydrogens to atoms CHARMM formally single-bonds (notably a sulfoxide S=O, e.g. DMSO, which obabel would otherwise "complete" with 2 extra H). Each atom's valence is now pinned in the generated MDL molblock (the V2000 ``vvv`` field = its bond-order sum), so ``obabel --gen3d`` returns exactly the RESI's atoms
- performance: the ``ring_check`` task's ``cutoff`` default is lowered from 10.0 Å to 4.0 Å (and likewise the ``RingChecker``/``ring_check`` library defaults). ``cutoff`` is only the candidate search radius -- each candidate is then gated by ``PSFRing.pierced_by``'s fixed 3.5 Å test -- so any value at or above 3.5 Å yields **identical** piercings while a smaller radius prunes far more candidate bond/ring pairs. On a ~200k-atom membrane the whole-system scan drops from ~38 s to ~13 s (and ~62 s to ~8 s on a periodic lipid patch). Verified result-identical at 3.5/4.0/10.0 on the known-piercing fixtures (1 and 2 piercings) and on the ex17 embedded membrane; a cutoff-invariance regression test is added

## [3.1.0] - 2026-07-07

- enhancement: ``modify-package pdb-repo add-entry`` now installs a whole **directory of entries** in one command, not just a single residue. ``make-pdb-collection --streamID <stream>`` writes a container directory with one ``<RESI>/`` subdirectory per residue; previously each had to be installed separately. Point ``add-entry`` at that container and it validates every entry up front, groups them by target collection (each residue's segtype default, or all into ``--collection NAME``), and repacks each affected tarball once. A single ``<RESI>/`` directory still works exactly as before -- the argument auto-detects which it is by the presence of ``info.yaml``. New ``ResourceManager.add_pdb_collection`` backs it (``add_pdb_entry`` remains the single-entry API); with ``--branch`` the commit message summarizes the batch
- new feature: the ``solvate`` task can now solvate with a **non-water solvent**. Its new ``solvent`` key (default ``TIP3``) selects the solvent species: ``TIP3``/``water`` uses VMD's built-in pre-equilibrated water box (unchanged), while any other RESI name is looked up as a ``kind: box`` entry in the ``solvent`` PDB collection -- pestifer checks out its pre-equilibrated box psf/pdb and passes them to VMD's ``solvate`` plugin as ``-spsf``/``-spdb`` with the box's exact equilibrated edge (``-ws``) and key atom (``-ks``), so that box tiles the cell instead of water. ``autoionize`` still runs afterward (``salt_con``/``cation``/``anion`` unchanged). A clear error is raised if the named solvent has no box entry (with the `make-pdb-collection solvent` command to build one) or is only a single-molecule entry. Validated by tiling a custom box onto a real protein solute via VMD
- enhancement: ``make-pdb-collection solvent`` can now build a box for a solvent that ships **no coordinates** -- neither a PDB-repository entry nor an internal-coordinate (IC) table -- by falling back to **Open Babel**. The RESI's atom+bond graph (atoms, elements, bond orders) is serialized to an MDL molblock and expanded to 3-D by ``obabel --gen3d``, then written out under the RESI's own atom names (RESI atom order is preserved, so the coordinate-to-name mapping is exact and needs no re-matching). This unlocks the common small-molecule CGenFF solvents (``MEOH``, ``ETOH``, ``DMSO``, …) that previously could not be built. The rough generated geometry only seeds the build; the minimize + NPT equilibration relaxes it. Requires ``obabel`` on PATH (already used by the ligand-parametrization path). Validated end to end: a ``MEOH`` box built via this path (edge ~20.4 A, density ~0.79 g/cc) solvates a protein solute cleanly (125 methanols placed)
- new feature: ``make-pdb-collection solvent`` builds a **pre-equilibrated periodic solvent box** for a single CHARMM residue -- the ``kind: box`` artifact VMD's ``solvate`` plugin needs to bulk-solvate with anything other than its built-in TIP3P water. `pestifer make-pdb-collection solvent --resname R --nmol 216 --density D` builds one molecule from the residue, packs `nmol` randomly-oriented copies into a cube sized for the target density, psfgens the box as a few segments of many residues each (like water in a membrane build), minimize + NPT-equilibrates it under PBC, then ships the *equilibrated* coordinates (NAMD `wrapAll` leaves every molecule whole and periodic-consistent, exactly what `solvate -ws` tiles cleanly) as an installable `<output-dir>/<RESN>/` entry (`info.yaml` with the equilibrated `box_edge`/`density`/`key_atom` + the box psf/pdb). Validated against VMD's 216-water reference: TIP3 -> ~18.6 A edge, ~1.0 g/cc, periodic-clean. Install it with `modify-package pdb-repo add-entry <dir> --collection solvent`. The existing flag-based conformer/stream sampler remains the default mode of `make-pdb-collection` (no verb), so `--streamID ...` invocations are unchanged
- enhancement: PDB-repository entries can now declare a ``kind`` in their ``info.yaml`` -- ``molecule`` (the default; a single-molecule conformer set, used by the grid membrane packer's chamber fill) or ``box`` (a pre-equilibrated periodic solvent box for VMD's ``solvate`` plugin, carrying ``box_edge``, ``key_atom``, and a psf/pdb pair). ``PDBCollection`` loads both kinds (a missing ``kind`` is treated as ``molecule``, so existing collections are unaffected), and ``PDBInput`` gains ``is_box``/``get_box_edge``/``get_key_atom``/``get_box_psf``/``get_box_pdb``
- change: the built-in PDB-repository ``water_ions`` collection is renamed to ``solvent`` (``water_ions.tgz`` -> ``solvent.tgz``; its 12 water/ion entries are unchanged). This is the first step toward a broader ``solvent`` collection that will also hold pre-equilibrated boxes of non-water solvents for VMD's ``solvate`` plugin. ``add_pdb_entry`` now routes ``ion``/``water`` residues to ``solvent`` by default. The CHARMM *stream* label ``water_ions`` (used to classify residues from ``toppar_water_ions.str``) is unchanged -- the pestifer collection and the CHARMM stream now intentionally differ
- change: example storage folders lost their ``ex`` prefix -- ``pestifer/resources/examples/ex01`` and ``docs/source/examples/ex01`` are now ``.../01``, ``.../02``, ... (the folder-name format is just the zero-padded ID). This also fixes a latent inconsistency where a newly-added example's generated ``.rst`` pointed its ``literalinclude`` at ``examples/NN/`` (no prefix) while the docs folder was created as ``exNN/``; both now agree
- change (**breaking**): the `modify-package` CLI is reorganized under a required **category** positional and a **verb**, replacing the flat `--example-*`/`--add-*`/`--regenerate-*` flags. The three categories are `example` (`add`/`update`/`delete`/`rename`/`author`), `pdb-repo` (`add-entry`), and `charmmff` (`add-residue`, `regenerate-segtypes`, `update-atomselect-macros`). So `pestifer modify-package --add-residue X.str` becomes `pestifer modify-package charmmff add-residue X.str`; `--example-action rename --example-id 19 --new-example-name N` becomes `example rename 19 N`; `--regenerate-segtypes` becomes `charmmff regenerate-segtypes`; and the new PDB-entry install is `pdb-repo add-entry DIR`. The shared `--branch` contribute flow (clean-tree check, branch, commit exactly the touched files) is available on the `pdb-repo` and `charmmff` verbs. `charmmff` is the home for future force-field-modification work
- new feature: `modify-package` can now contribute a **PDB-repository entry** -- coordinates for a residue so `make_membrane_system` can place it -- completing the "contribute buildable content for review" workflow alongside `charmmff add-residue`. `pestifer modify-package pdb-repo add-entry DIR` takes an entry directory produced by `make-pdb-collection` (named after the residue, containing `info.yaml` plus its conformer PDBs), validates it (info.yaml parses; every referenced conformer PDB exists), and installs it into the built-in PDB repository. A repository is a set of collections, each a `<stream>.tgz` tarball holding one `<RESI>/` subdirectory per residue; the entry is added under `<collection>/<RESI>/` (existing entries preserved), creating the collection tarball if it does not exist. `--collection NAME` picks the target (default: the residue's segtype, with `ion`/`water` -> `solvent`); a resname already present is refused unless `--force`. Combined with `--branch NAME` it commits exactly the changed collection tarball on a fresh branch, like the residue-contribution flow. `add_pdb_entry` also accepts `kind: box` entries (a solvent box's psf/pdb pair rather than conformer PDBs), installing them into a collection the same way

## [3.0.1] - 2026-07-07

- bugfix: VMD steps silently did nothing on hosts whose `vmd` launcher wraps the binary in `rlwrap` (common in site VMD 2.x module builds). Pestifer launched every external command with `start_new_session=True` (so a NAMD `charmrun`->`namd3`->PE tree can be torn down with one `killpg` on interrupt), but a new session strips the controlling terminal, and an rlwrap-wrapped VMD then exits ~immediately (returncode 0, no output) without ever running its Tcl script. The visible symptom was `desolvate` failing its output guard ("VMD returncode 0 ... expected outputs were not produced") even though `vmd` ran fine by hand. VMD is now launched **in pestifer's own session** (`new_session=False`) and signaled by pid rather than process group, and its stdin is `/dev/null` so such a launcher never selects rlwrap in the first place; NAMD/charmrun keep their own session and `killpg` teardown. The `desolvate` guard message now distinguishes "VMD did not launch" (rc != 0) from "VMD ran but produced no output" (rc 0), instead of always blaming a missing `vmd`

## [3.0.0] - 2026-07-06

- change (**breaking**): the **packmol** membrane packer has been removed; the grid packer is now the only way to build a bilayer. `make_membrane_system` placed lipids either with packmol (small "patches" tiled into a "quilt") or on a deterministic grid; the grid packer is orders of magnitude faster, needs no external binary, and has been the default in every shipped membrane example, so packmol is retired. The `bilayer.packer` option now accepts only `grid` (a config with `packer: packmol` gets a clear validation error), and the packmol-only keys (`nloop`, `nloop_all`, `tolerance`) and the `paths.packmol` executable setting are gone. **packmol is no longer a runtime dependency** -- pestifer no longer checks for or requires a `packmol` executable at startup. Removed with it: the packmol scripter and log parser, the `bilayer_quilt.tcl` tiling script, and the packmol-specific artifact types. The asymmetric-membrane build still constructs two symmetric calibration "patches" to measure per-leaflet preferred areas, and the `relaxation_protocols` `patch`/`quilt` keys are retained (the grid path uses both). Existing configs that set `packer: grid` (including examples ex16/ex17) are unaffected
- bugfix: the `desolvate` task no longer fails silently and misleadingly when its VMD index/PSF-generation step does not run. That step writes the atom-index file (`dry.idx`) and dry PSF that `catdcd` then consumes; if VMD does not actually execute -- e.g. `vmd` is not launchable on the invoking shell's PATH, in which case it can exit 0 having produced nothing -- the index was never written, `catdcd` failed downstream with a confusing `Error opening index file` / `Error reading index file`, and the task still reported `result: 0`. `DesolvateTask` now verifies the index and PSF were actually produced (and the VMD return code) before invoking `catdcd`, and aborts with a clear message pointing at the real culprit (VMD not running) plus a one-line check to reproduce it; it also now checks `catdcd`'s own return code
- bugfix: standalone post-processing subcommands (`desolvate`, `mdplot`) are no longer hard-gated on the full build toolchain. `Config._set_shell_commands` unconditionally required `charmrun`, `namd3`, `vmd`, and `catdcd` (before their removal in this release, also `packmol`) to all resolve on PATH and aborted otherwise, so e.g. `pestifer desolvate` (which needs only `vmd`/`catdcd`) refused to run unless the full NAMD toolchain was also loaded. Missing commands are now fatal only for build runs (a user config file was supplied, `verify_access=True`); otherwise they are recorded and warned about, and the task that actually uses a missing command fails loudly itself when it needs it

## [2.9.0] - 2026-07-06

- bugfix: the pre-execution pipeline validation (2.8.0) fired on a controller configured with **no user tasks** -- `configure()` auto-appends a `terminate` task, and validating that lone terminate reported "task 00 'terminate' requires an existing molecular system." This path is taken when a controller is configured empty and its task list is supplied later via `reconfigure_tasks` (and would also have rejected a genuinely empty build with a confusing message). Validation now runs only when the user actually supplied tasks; an empty task list is a placeholder/no-op, not a pipeline
- change: residue **segtype classification is now derived from the CHARMM force field** instead of being hand-maintained. Previously `core/labels.py` hard-coded ~660 residue names across segtype lists (the lipid and glycan lists alone were ~360 names that had to be extended by hand whenever CHARMM added residues). These are now derived from the topology/stream file each residue is defined in (`*_lipid_*` -> lipid, `*_carb_*` -> glycan, `*_prot_*` -> protein, `*_na_*` -> nucleicacid, `*_cgenff*` -> ligand, `water_ions` -> water/ion), persisted to a generated resource (`resources/labels/derived_segtypes.json`) and merged into `Labels.segtype_of_resname` at import. `labels.py` now holds only the ~120 *curated* names that cannot be derived (PDB/IUPAC sugar codes such as `NAG`/`MAN`/`GLC` which have no CHARMM RESI, cross-category exceptions such as `HEME`/`NAD`/`ADP`, the standard protein residues used to expand psfgen pdbalias wildcards, and the water models used to split water from ions). The change is classification-preserving -- every residue that was classified before keeps its exact segtype (verified against the pre-change table) -- while **4,000+ additional force-field residues are now classified**, so a lipid or glycan added to CHARMM's stream files is recognized without editing pestifer. Regenerate the derived resource after a force-field update with `pestifer modify-package --regenerate-segtypes`. (The VMD atomselect macros still build from the curated lists; broadening their scope is deferred.)
- change: `ring_check` now removes its rotation-resolution **scratch files** after a piercing is cleared. The aromatic side-chain and glycan-pendant resolution paths write a trial PDB per candidate angle (each a full-system coordinate file) plus the VMD generator script/log; only the winning pose is copied to the registered `ring_check-NN.pdb` artifact, so the rest were pure litter left in the run directory outside the artifacts tarball (a real membrane build left ~17 files / ~67 MB). They are now deleted once the accepted pose is registered, and retained only on the failure path (for debugging a piercing that no rotation could clear)
- new feature: `modify-package` can now contribute a **built-in custom residue** for review. `pestifer modify-package --add-residue FILE` validates a CHARMM topology file (`.str`/`.rtf`/`.top`), extracts the `RESI` names it defines, refuses names that already exist in the force field (unless `--force`), copies the file into the active force field's `custom/` directory, registers each `RESI` name under a segtype in `core/labels.py` (default `ligand`, overridable with `--segtype`), and clears the resource cache so the new residue is picked up on the next run. The `--branch NAME` option folds the git workflow into the command: it requires a clean working tree, creates a new branch off the current HEAD, makes the change, and commits **exactly** the files it touched (the copied custom file and, if a segtype was registered, `labels.py`) -- never anything else in the tree -- then prints the `git push`/`gh pr create` commands to open a pull request. Pushing and PR creation are left to the contributor. The segtype registration edits the `_segtypes` literal in `labels.py` in place via an `ast`-located, format-preserving source insertion that is re-parsed before write-back
- bugfix: re-running `make_membrane_system`'s embed step in a directory that still holds intermediates from a prior run silently produced an **unsolvated** system -- no water above (or below) the bilayer around the protein. In the no-salt path the water slab is promoted to `*_solution_{upper,lower}.{psf,pdb}` with a plain `file copy`, which throws `file already exists` when a stale copy from the previous run is present; that error fired *before* the slab was appended to `addl_solution`, and because `addl_solution` was never initialized (the list was created under the unused name `addl_water`), the append loop then died on `can't read "addl_solution": no such variable` and was skipped entirely -- yet the embed still wrote its output and reported success. The missing water leaves a tall vacuum column around the protein that passes minimize and NVT (both fixed-box) and then collapses on the first NPT step with `Periodic cell has become too small for original patch grid`. `bilayer_embed.tcl` now initializes `addl_solution`, uses `file copy -force` for all solution and final-output copies so a re-run does not collide with stale intermediates, and hard-fails (rather than shipping an unsolvated system) if solvate generated slabs but none were appended

## [2.8.0] - 2026-07-06

- new feature: the task pipeline is now validated **before** execution, so malformed task hand-offs fail in milliseconds with an explanation instead of crashing partway through a build (or silently discarding work). Each task declares a contract -- the pipeline "currencies" it requires and provides (`base_coordinates` from `fetch`; the built `state` fileset; the in-memory `base_molecule` object that only the `psfgen` family produces; `md_output`) -- and a pre-execution check walks the task list once, reporting: a task that needs a currency no earlier task supplies (e.g. `md`/`solvate` as the first task with no system; `continuation` → `psfgen`, which would rebuild from scratch and discard the continued system; `cleave`/`ligate`/`pdb2pqr` after a `continuation`, which provides files but not the `base_molecule` they need; `mdplot` with no preceding `md`; a dangling `fetch`), an origin task that would discard an already-built system, and anything scheduled after the terminal `terminate`. Double-solvation is a warning rather than an error. `desolvate` is recognized as a standalone post-processing utility (it reads/writes explicit file paths, not pipeline currencies) and is rejected with an explanation if placed in a build `tasks:` list, while `pestifer desolvate` and `pestifer mdplot` continue to run standalone (they skip the check)

## [2.7.2] - 2026-07-02

- new feature: `ring_check` now **resolves** pierced rings rather than only reporting them, choosing the action by the pierced ring's segtype. `ring_check.segtypes` has **no default** and must be set explicitly (e.g. `[lipid, glycan, protein]`; an empty list checks nothing), and the cycle search is length-bounded (`max_ring_size`, default 7) so the giant chordless cycle a disulfide makes through the backbone is skipped while real 5/6-membered rings are found. An **aromatic** protein side-chain ring (His/Phe/Tyr/Trp) is swung off the piercing bond by rotating the side chain (chi2 then chi1, via the `PestiferCRot` `SCrot_chi2`/`SCrot_chi1` procs, whose `SCrot_chi2` coordinate-parsing bug is fixed here). A **rigid** ring that cannot itself move -- a **proline** side-chain ring (fused to the backbone) or a **glycan** ring -- when speared by a glycan bond is cleared by rotating the glycan sub-branch that carries the piercing bond about an upstream rotatable bond (a *bridge* of the glycan graph), swinging the bond out as a rigid body; this rotation is done in memory (numpy) and only the single winning pose is written. Among the rotations that un-thread a ring, the one introducing the fewest new close heavy-atom contacts is kept, so a pierce fix does not trade a thread for a steric jam. A **lipid** ring is resolved by deleting a segment (`delete`: `piercee`, `piercer`, `both`, or `none`); `both` removes the piercee and the piercer when it too is a lipid -- robust for inter-threaded lipids, which strain both molecules. The topology (PSF parse, bond list, ring cycles) is parsed once into a reusable `RingChecker` and only coordinates are re-checked per candidate. (The membrane example ex17 now sets `segtypes: [lipid, glycan]` explicitly, where it previously relied on the removed default.)
- enhancement: the whole-system pierced-ring scan is vectorized -- bond and ring coordinates come from a single positional coordinate array (in PSF atom order) bucketed into a link cell, instead of a per-element pandas `.loc` ingest -- so a several-hundred-thousand-atom membrane scans in ~30 s instead of minutes. The positional read also fixes a `KeyError` on systems over 100k atoms, whose PDB atom serials wrap past 99999 and so no longer match the PSF's 1..N serials
- bugfix: a grid-packed membrane's condensing relaxation minimize can pull a lipid tail through a sterol ring; the threading is topological (further minimization does not undo it) and RATTLE-fails the first dynamics step. `make_membrane_system` now runs a lipid `ring_check` -- followed by a short minimize to relax the freed lipid and the deletion void -- between the opening minimize and the first dynamics stage of a grid build's relaxation, so the pierce is removed before dynamics. Only grid builds are guarded; packmol places lipids with a hard overlap tolerance
- bugfix: after embedding, the water slabs added above/below the protein to fill the z-expanded box were solvated from the box origin (`{0 0}` to `{Lx Ly}` laterally), which assumed the membrane occupied the +x/+y quadrant starting at the origin. After NAMD relaxation the membrane is wrapped/re-centered away from the origin, so the added water came out offset by tens of A in x and y. The slabs are now centered on the membrane (`mid +/- L/2`) so they register with it
- bugfix: the grid packer builds the lipid and chamber-solvent lattices independently, so a solvent molecule could land almost coincident with a lipid atom (observed: a TIP3 0.39 A from a PSM choline methyl). Such a near-coincidence makes VMD infer spurious bonds (`Exceeded maximum number of bonds`) and psfgen then fails to read that lipid residue's coordinates and silently guesses them (a mis-placed hydrogen). `write_grid_pdb` now indexes the placed lipid atoms in a cKDTree and drops any chamber-solvent molecule within `clash_cutoff` (1.2 A) of a lipid -- a handful of waters, negligible for hydration, while the relaxation still resolves the milder overlaps
- change: secondary-structure restraints (`md.ssrestraints`) are now opt-in via `ssrestraints.enabled` (default off) instead of being built for every md stage. Every attribute of the `ssrestraints` schema dict was defaulted, so ycleptic auto-populated it and the md task ran the VMD ssrestraints plugin on *every* stage -- including bare-membrane bilayer-relaxation and calibration-patch stages, which have no secondary structure (each spun up VMD to select nothing and wrote an empty extraBonds file). Protein/nucleic stages that should hold their secondary structure must now set `ssrestraints: {enabled: true}`

## [2.7.1] - 2026-07-01

- bugfix: the embed step (and the bare-membrane output path) used the *pre-relaxation* bilayer coordinates instead of the relaxed ones, silently discarding the quilt relaxation. `equilibrate_bilayer` relaxes the quilt in a subcontroller and rekeys the relaxed state to `<name>_state` before merging it back, but `import_artifacts` kept the parent pipeline's existing same-key head entry (the pre-relaxation psfgen state registered by `do_psfgen`) and buried the incoming relaxed one, so `embed_protein` read un-relaxed coordinates. This affected both packers; the grid packer made it obvious because its pre-relaxation lattice is loose. `import_artifacts` now lets an imported head artifact supersede a stale same-key entry (the superseded entry is moved to history), so the relaxed bilayer becomes current
- bugfix: grid lipid placement (`write_grid_pdb`) filled a fixed nx×ny lattice row-major, leaving the trailing `nx·ny − n` positions empty as a lipid-free strip at the box edge; lipids are now distributed evenly across the rows and spaced across the full box width so each leaflet is gap-free
- bugfix: a non-embedding grid membrane build raised `IndexError` when computing `patch_nlipids × npatch`, because `npatch` (a schema list with no default) materializes as an empty list; it now defaults to `[1, 1]` when `npatch` is unset, empty, or malformed
- docs: the `make_membrane_system` documentation now covers the grid packer alongside packmol -- a packer-agnostic intro, the `packer` option, separate grid/packmol task-flow descriptions and diagrams, and updated relaxation/SAPL notes. Membrane-building mentions elsewhere (introduction, installation, `make-pdb-collection`, the `make_membrane_system` schema description, and the README) no longer imply packmol is required, and a note clarifies that a grid build's `*-quilt` artifacts are the full-size membrane (a legacy term), not patches

## [2.7.0] - 2026-06-25

- bugfix: the `Bilayer` constructor filled per-species lipid/solvent counts (`patn`) into the caller's `composition_dict` **in place**. Because `make_membrane_system` builds a small calibration/initial patch (default 100 lipids/leaflet) from that dict before deep-copying it for the full gridded membrane, the patch's counts got baked in and reused: the full membrane was populated with only ~100 lipids and hydrated for only ~100 lipids regardless of its actual (often much larger) lipid count. Under-population + under-hydration (water-per-lipid far below the requested ratio) made the membrane condense far past its tensionless area into a gel-like collapse during relaxation, and left the embedded system grossly over-volume so production could not reach proper density. The constructor now works on its own copy, so each membrane computes counts from its own `leaflet_nlipids`. The symptom was masked at small box sizes (where ~100 lipids happened to match the protein-footprint box) and only surfaced once larger boxes needed more lipids
- new feature: `NPgT` is now the preferred name for the constant-x:y-ratio barostat ensemble (`useFlexibleCell` + `useConstantRatio`: x and y scale together while the cross-sectional area relaxes, z independent), useful for membranes. `NPgT`/`NPGT`/`npgt` are equivalent. Pestifer historically called this `NPAT`, but that label has never set `useConstantArea` (true constant area) -- only constant ratio -- so the "A = area" in the acronym was misleading. `NPAT` is still accepted and behaves identically to `NPgT`, but now emits a one-time deprecation warning: a future version will reserve `NPAT` for genuine constant-area simulations. Affects the `md` task ensemble and every membrane-relaxation protocol stage
- new feature: the `make_membrane_system` bilayer builder gains a grid packer (`bilayer.packer: grid`, alongside the default `packmol`) that places lipids and solvent on a lattice deterministically instead of running packmol's constrained packing — far faster, and with no packmol dependency (`paths.packmol` is not needed for grid builds). It builds the full membrane directly, sizing the box to the embedded protein's xy footprint plus the `embed.xydist` margin when embedding (or to `patch_nlipids` × `npatch` otherwise), so the patch→quilt tiling is bypassed entirely for grid builds
- new feature: asymmetric membranes (compositionally different leaflets) are built stress-free by the grid packer without leaflet extraction. Two symmetric "calibration" patches — one per leaflet composition — are relaxed to measure each leaflet's preferred area-per-lipid (APL), and the full membrane is then gridded at per-leaflet counts in the stress-free ratio `n_upper/n_lower = apl_lower/apl_upper` (sized to the protein footprint when embedding), giving zero differential stress by construction. The patch→quilt path now serves only the packmol packer
- new feature: an optional `diagnose_differential_stress` block on `make_membrane_system` runs a short tensionless (NPγT) pressure-profile pass on an assembled asymmetric membrane and reports the per-leaflet surface tensions, their difference (the differential stress Δγ), and the lipid-count move that would null Δγ (estimated from an assumed area-compressibility modulus). Diagnostic only — the membrane is not rebuilt; requires a CPU (non-CUDA) NAMD build
- enhancement: asymmetric calibration warns when a calibration patch has not equilibrated (its lateral area is still drifting over the final half of relaxation), since an under-relaxed cell yields an untrustworthy preferred APL — the classic symptom is two differently-composed leaflets reporting near-equal APLs because a cholesterol-rich leaflet failed to condense
- enhancement: long NPT/NPgT stages of a gridded-membrane relaxation are automatically split into a ramp of restarts, so the cell can condense (a gridded membrane is built at the extended lipid length and shrinks tens of percent under a barostat) without any single run shrinking it past NAMD's startup spatial-decomposition patch grid (`Periodic cell has become too small for original patch grid`). Restarting is NAMD's recommended remedy -- a fresh run rebuilds the patch grid for the smaller cell -- and avoids the slower `margin` workaround. Sub-runs are kept at or above the Langevin-piston period so the barostat actually moves the cell
- bugfix: grid placement clamps `half_mid_zgap` to a safe minimum. Deterministic placement puts opposing-leaflet tail atoms at `midplane_z ± half_mid_zgap`, so a `half_mid_zgap: 0` (harmless for packmol, which resolves overlaps with its tolerance) made them coincide and the relaxation minimize diverged with a segmentation fault
- bugfix: the asymmetric grid build sizes its loose starting box from the calibrated APL (`1.15 × max(apl_upper, apl_lower)`) rather than `max(SAPL, 1.15 × apl)`; with a condensed (cholesterol-rich) leaflet the nominal `SAPL` badly overestimates the real APL, gridding the box so loose it had to shrink ~25% and overflowed the patch grid
- bugfix: the `mdplot` task no longer aborts when a pressure-profile block has a single profile row — its NaN standard deviation (ddof=1) made the min/max-with-σ x-range NaN and `set_xlim` raised "Axis limits cannot be NaN or Inf"; the range is now computed nan-safely and explicit limits are set only when finite
- change: examples 16 (symmetric DMPC) and 17 (asymmetric model viral bilayer) now build their membranes with the grid packer (`packer: grid`) instead of packmol
## [2.6.3] - 2026-06-16

- bugfix: interrupting a build (Ctrl-C, or `kill`/scheduler signal) now shuts down cleanly instead of orphaning the running NAMD job. External commands are launched via `subprocess.Popen(..., shell=True)`, so the tracked pid was the shell/`charmrun` wrapper; a `SIGTERM` to it never reached `namd3`, leaving an N-PE NAMD process pinning every core after pestifer exited (the old `atexit`-based `kill_child` also never fired on a signal, only on a clean exit, and tracked a single module-global pid). Each external command now runs in its own session (`start_new_session=True`) so its whole process tree (`charmrun` → `namd3` → PEs) can be removed with one `killpg`; pestifer installs `SIGINT`/`SIGTERM` handlers that terminate every live child group (SIGTERM, then SIGKILL stragglers), log a clear `pestifer run INTERRUPTED by <signal>` message naming the working directory, and exit with status `128+signum`. Handler installation is a no-op off the main thread and under pytest

## [2.6.2] - 2026-06-16

- bugfix: the `merge` task (and the `continuation` task) no longer abort when an input PSF was built outside pestifer and records a topology stream file pestifer does not ship. VMD writes `REMARKS topology .../toppar_water_ions_namd.str` (its NAMD-specific variant of the standard `toppar_water_ions.str`, which pestifer *does* ship and which covers the same residues); `merge` harvested that name from the PSF REMARKS and emitted a psfgen `topology toppar_water_ions_namd.str` line for it, but `copy_charmmfile_local` only warns and writes nothing for an unknown file, so psfgen died with `ERROR: Unable to open topology file ...` → `MOLECULE DESTROYED BY FATAL ERROR`. Both tasks now keep only the stream files that actually resolve locally (after the copy attempt), drop the rest with a clear per-file warning naming the file and the source PSF, and `merge` additionally strips the dropped file's carried-over `REMARKS topology` line from the merged PSF so a downstream task reading it cannot re-trigger the same failure. Dropping is safe because the merge is a pure `readpsf` of fully-built input PSFs (and the standard water/ions stream is still loaded)

## [2.6.1] - 2026-06-11

- bugfix: large membrane-embedding builds aborted at the end of the `make_membrane_system` embed step with a VMD segmentation fault (signal 11), even though the embedded, solvated, and neutralized system had already been written correctly. The clash-removal steps in `bilayer_embed.tcl` built `atomselect "index <list>"` selections directly from the per-atom index lists returned by `measure contacts`; on a large solvated system those lists run to ~10^5 atoms, and VMD's atom-selection parse tree has one node per index whose destructor (`atomparser_node::~atomparser_node`) is recursive — so freeing such a selection recurses ~10^5 frames deep and overflows the stack during interpreter teardown at VMD exit (confirmed via core dump: `VMDApp::~VMDApp` → `DeleteInterpProc` → `AtomSel::~AtomSel` → `ParseTree::~ParseTree` → hundreds of thousands of nested `~atomparser_node` frames). The crash was purely in at-exit teardown but returned signal 11 and aborted the whole build, and it reproduced on VMD 1.9.3, 2.0.0, and 2.0.1a1 (it is a long-standing VMD bug, not version-specific). The clashing atoms are now mapped to `(segname, resid)` by reading directly from the parent whole-molecule (`all`) selection's attribute lists, so no large `index` selection — and no large parse tree — is ever constructed
- bugfix: the same clash-removal loops invoked `delatom $seg $resid` (which deletes an entire residue) once per *contacting atom* rather than once per residue, so a lipid or water with several atoms near the protein was deleted on the first call and then re-deleted on every subsequent atom — emitting tens of thousands of `no residue ... of segment ...` messages on a large embed. Each clashing residue is now deduplicated and deleted exactly once
- enhancement: `embed_protein` no longer appends a spurious trailing `regenerate angles dihedrals` to the embed psfgen script. It called `writescript(regenerate=True, writepsf=False, writepdb=False)`, which emitted a `regenerate` with no following `writepsf`/`writepdb` to consume it — pure wasted work regenerating angles and dihedrals for the ~2M-atom post-autoionize context. The `bilayer_embed` script already regenerates before each structure it writes, so `embed_protein` now uses `regenerate=False` (matching the sibling quilt call site)

## [2.6.0] - 2026-06-10

- new feature: NAMD runs under SLURM can now span multiple nodes, controlled by a new `namd.cpu-parallel-launcher` option (`auto` | `numactl` | `srun` | `mpirun` | `charmrun`, default `auto`). Previously the SLURM CPU launch path always used `numactl --interleave=all namd3 +p N`, which runs a single node-local process and cannot use more than one node; on a multi-node allocation with an MPI NAMD build this failed immediately with the Charm++ "to use multiple processors ... charmrun +pN namd3 / build the mpi-linux-x86_64-smp version" error. `auto` keeps the single-node `numactl` launch for one-node allocations and switches to `srun namd3` when `SLURM_NNODES > 1`; `mpirun` launches via the MPI runtime's own process manager (`mpirun -np N namd3`, e.g. Intel MPI Hydra); `charmrun` uses charmrun's `++mpiexec` launcher for Charm++ net-layer builds. A companion `namd.srun-mpi-type` option emits `srun --mpi=<plugin>` (e.g. `pmi2` / `pmix`) for MPI runtimes — notably Intel MPI — whose ranks would otherwise each initialize as an independent 1-PE job and collide on the output DCD with `Couldn't open DCD file ...: File exists`
- bugfix: multi-node NAMD launches no longer stage parameter files to node-local scratch. `_stage_params_to_local_scratch` copied the parameter files to `$TMPDIR` and rewrote the NAMD config's `parameters` line to that absolute path, which is only valid when every PE runs on one node; under a multi-node `srun`/`mpirun` launch the files existed only on the launch node, so ranks on the other nodes died with `UNABLE TO OPEN CHARMM PARAMETER FILE`. Staging is now gated on a single-node launch, and multi-node runs read the parameter files from the shared filesystem (as they already do for `structure`/`coordinates`)
- bugfix: the NAMD log parser no longer aborts the whole build on a single malformed line. Live NAMD stdout from a multi-rank run (e.g. `srun`/`mpirun` across nodes) can interleave output from different PEs, producing a garbled `ENERGY:` record with a token like `0LINE`; `process_energy_line` then raised `ValueError` on `float()` and killed the run even though the simulation itself was fine. `process_line` now catches `ValueError`/`IndexError`/`KeyError` from a single line processor, logs it at debug level, and continues without appending a partial record
- bugfix: subcontroller tasks now inherit the parent run's NAMD settings. `Config.taskless_subconfig` built the subcontroller configuration from schema defaults only, so MD runs launched inside `make_membrane_system` (which uses a subcontroller for the membrane-relaxation stages) ignored user overrides such as `namd.cpu-parallel-launcher` and reverted to the default — reintroducing the multi-node launch failure for the relaxation stages even when the top-level tasks were correctly configured. The subconfig now inherits the progenitor's user settings (NAMD launcher, paths, force field, psfgen options) and resets only the task list
- bugfix: `terminate`/`package` builds that run no MD step (e.g. `continuation` → `psfgen` → `terminate`) produced a consolidated minimal `.prm` with an empty NONBONDED section (`nonbonded=0`) and only partial bonded terms, so the packaged production run failed with `DIDN'T FIND vdW PARAMETER FOR ATOM TYPE ...` (the first such type encountered, e.g. `NH3`). The standard CHARMM `.prm` parameter files are staged by NAMD/MD tasks, so a build with no MD step had registered only topology `.str` stream files for consolidation; `generate_minimal_params` now always stages and includes the standard parameter set (via the NAMD scripter's `fetch_standard_charmm_parameters`) before consolidating, and warns loudly if a consolidated file still has no nonbonded parameters
- bugfix: corrected a Tcl variable-dereference typo (`set $box_min_z` / `set $box_max_z` instead of `set box_min_z` / `set box_max_z`) in the bilayer-embed water-gap padding, which created a garbage-named variable and left the box bound unchanged when the water gap above or below the protein was under 3 Å
- dependency: dropped the `gputil` runtime dependency (and the `setuptools` shim it required to provide `distutils` on Python 3.12+). Local NVIDIA GPU enumeration now queries `nvidia-smi --query-gpu=index` directly via a small helper, with the same behavior as before (enumerates physical GPUs, returns an empty list when `nvidia-smi` is absent or fails)

## [2.5.2] - 2026-05-28

- bugfix: a plain `pip install pestifer` (without the optional `ligand-paramgen` extra) crashed on every CLI invocation with `ModuleNotFoundError: No module named 'dimorphite_dl'`; the optional `rdkit` and `dimorphite_dl` dependencies were imported at module load time in `charmmff.ligand_paramgen.protonation` and `charmmff.ligand_paramgen.mol2_writer`, and since `subcommands/__init__.py` eagerly imports the `make-ligand-mol2` subcommand at startup, those imports ran for every command. The imports are now lazy (deferred into the functions that use them) and raise a clear `pip install pestifer[ligand-paramgen]` hint when actually needed
- bugfix: the `dependencies` array in `pyproject.toml` was missing its closing `]`

## [2.5.1] - 2026-05-28

- bugfix: `PSFContents.remove_ignored_residues` crashed with `TypeError: 'PSFAtom' object is not subscriptable` because it wrongly indexed the result of `BaseObjList.get()` with `[0]`; `get()` returns the matched object directly when there is exactly one match (which is always the case for an atom-by-serial lookup), so the trailing `[0]` was indexing into the atom itself
- bugfix: `AsymmetricUnit` extended `ignored_residues` twice with the same `more_ignored_residues` list (a copy-paste duplicate around the intervening `grafts.assign_residues` call), causing `PSFContents.remove_ignored_residues` to be called for the same residues twice and crash on the second pass when `self.atoms.get()` returned `None` (atoms already removed) and `self.atoms.remove(None)` invoked `PSFAtom.__eq__` against `None.serial`
- bugfix: continuation+psfgen builds that inserted protein residues at glycosylation sites (e.g. extending the HIV-1 gp41 ectodomain from MPER residue 685 through TM residue 710 on a pre-built glycoprotein) emitted glycan-link patches against the protein segment (`patch 14ba A:1222 A:1223`) instead of the glycan segment (`patch 14ba AG01:1222 AG01:1223`), causing psfgen to abort with `no residue X of segment A` / `MOLECULE DESTROYED BY FATAL ERROR`; two related issues were fixed: (a) `Link._from_psflinkpatch` populated `chainID1`/`chainID2` from the psfgen segname but left `segname1`/`segname2` unset, so when `chainID1` was later remapped from the segname (`AG01`) to the biological chainID (`A`) so the lookup against PDB-derived residues could succeed, the writer at `PsfgenScripter.write_link` (which prefers `L.segname1` over `L.residue1.chainID`) had nothing but the biological chainID to fall back on; (b) `ResidueList.renumber`, when renaming a glycan residue whose resid collided with a newly inserted protein residue (glycan A:686 → A:1324 after inserting protein A:686), only propagated the resid update to links keyed by biological chainID and skipped links from PSF REMARKS whose `chainID1`/`chainID2` carried the psfgen segname, so those stale links became orphans and the orphan cleanup at `LinkList.assign_residues` removed the entire glycan tree from the residue list — leaving segment AG01 with only the renumbered root and stripping the patches' target residues; `renumber` now also keys its `mapper_by_chain` by segname when segname differs from chainID

## [2.5.0] - 2026-05-18

- new feature: `pestifer make-ligand-mol2` subcommand that, given a 4-letter RCSB PDB code or a path to a local `.pdb` file, generates a CGenFF-ready Tripos mol2 file for every HETATM ligand whose resname is not already defined in the loaded CHARMM topologies; for each unknown ligand it fetches the SMILES from the RCSB Chemical Component Dictionary (preferring the OpenEye canonical-stereo descriptor), builds an RDKit `Mol` whose heavy atoms are in PDB order and carry the original PDB atom names, protonates at the target pH (default 7.4) via Dimorphite-DL, places explicit Hs in 3D, and writes the mol2 via Open Babel with `--title <RESNAME>` so the CGenFF binary/web tool registers the resulting `RESI` under the correct residue name (rather than the leaked tempfile path); per-ligand failures are surfaced in a summary table without aborting the run, and `--smiles RESNAME=...` overrides the RCSB lookup when a custom tautomer or charge state is needed
- new feature: new `pestifer.charmmff.ligand_paramgen` subpackage exposes the building blocks reused by `make-ligand-mol2` as a library: `fetch_ligand_smiles(comp_id)` (RCSB lookup with LRU caching and HTTP-error → `RCSBLookupError` mapping), `protonate_ligand(residue, smiles, ph)` (RDKit + Dimorphite-DL pipeline that preserves PDB heavy-atom order and stamps each heavy atom with the original PDB name on the `_pdb_name` property), `write_mol2(mol, resname, outpath)` (PDB-block-via-obabel writer that propagates `_pdb_name` and assigns sequential `H1`, `H2`, ... names to hydrogens), and `generate_ligand_mol2s(parsed, charmmff_content, outdir, ph, smiles_overrides)` (the orchestrator the subcommand calls, which the build workflow will also call for unknown ligands)
- new feature: well-known per-user toppar cache at `~/.pestifer/toppar/`; when this directory exists, pestifer automatically prepends it to `charmmff.user_custom.searchpath` so `.str` files dropped there are picked up by every build without requiring a `user_custom` block in YAML; any project-local searchpath entry is loaded after the cache and overrides it (last-loaded wins) with a warning naming both files; designed as the default destination for `.str` files produced by running `make-ligand-mol2` output through the CGenFF binary or web tool
- enhancement: `charmmff.user_custom.searchpath` entries now expand `~/` (the user's home directory) and shell environment-variable references like `$HOME/foo` or `${SCRATCH}/charmm` before the directory check; previously those were taken literally and silently failed
- behavior change: `CHARMMFFContent.add_custom_directory` no longer asserts on duplicate-basename clashes between searchpath entries; instead, the later-loaded file wins and a warning is emitted naming both paths (the assertion previously crashed pestifer when, e.g., the user's project happened to contain a `.str` file with the same basename as one in the cache); same-resname clashes across differently-named files in the searchpath also emit a warning and follow last-loaded-wins semantics, so a project-local override of a cached parameterization is surfaced without silently merging
- new optional dep: `[project.optional-dependencies]` table in `pyproject.toml` declares an `ligand-paramgen` extra (`pip install 'pestifer[ligand-paramgen]'`) that pulls in `rdkit` and `dimorphite_dl`; the `obabel` binary remains a runtime-PATH dependency, documented alongside the new subcommand
- test: `tests/unit/test_charmmff/test_charmmffcontent.py::TestCharmmffContent.setUpClass` was picking the CHARMMFF version directory by `mtime` rather than by parsed month-year (the production `ResourceManager.charmmff_version_dirs()` orders semantically); when the `feb26` and `jul24` directories were extracted out of release order, `mtime` selected `jul24` and the test's hard-coded counts (calibrated against `feb26`) failed by ~1400 residues; setUpClass now delegates to `ResourceManager.charmmff_version_dirs()` so the test always exercises the same CHARMMFF release that production code uses
- test: new `tests/unit/test_charmmff/ligand_paramgen/` package covers the ligand-paramgen building blocks (29 tests): protonation invariants on five SMILES across a pH=7.4 charge range, RCSB descriptor-priority selection under `unittest.mock` of `requests.get`, mol2 title/atom-name round-trip via the real `obabel` binary, and orchestrator-level altloc dedup and known-resname filtering using real `Hetatm` objects; `test_charmmffcontent.py` gains three `TestAddCustomDirectory` cases (tilde expansion, duplicate-basename warning + last-wins, duplicate-RESI warning + last-wins) that build a stub `CHARMMFFContent` via `__new__` to avoid paying the full tarball-load cost
- doc: new `docs/source/subs/make-ligand-mol2.rst` covers the subcommand with examples (PDB code, local file, SMILES override, pH override) and outlines the CGenFF round-trip workflow (web tool → `~/.pestifer/toppar/`); linked from the subcommand toctree in `usage.rst`; the `charmmff.user_custom` schema docstring gains "Per-user toppar cache" and tilde-expansion sections

## [2.4.13] - 2026-05-14

- behavior change: `make_membrane_system` no longer silently turns on NAMD's on-the-fly `pressureProfile` calculation during NPT/NPAT bilayer-relaxation stages on CPU NAMD; pressure profile calculation is now opt-in via a new top-level `make_membrane_system.compute_pressure_profile` block (`enabled: false` default, plus `slabs` / `freq` knobs that map to `pressureProfileSlabs` / `pressureProfileFreq`); previously, any NPT/NPAT stage whose `other_parameters` did not explicitly mention `pressureProfile` had it auto-injected as `on` on CPU builds, surprising GPU users (the `processor-type` guard correctly skipped injection on GPU but the opt-out default was hostile and unintuitive); the `mdplot` task's `profiles=['pressure']` list, which was unconditionally requested by `equilibrate_bilayer`, is now likewise gated on the new flag so no spurious pressure-profile plotting is attempted when the data is not produced
- new feature: `make_membrane_system` now hard-errors with `PestiferBuildError` if `compute_pressure_profile.enabled: true` is set under GPU NAMD (`processor-type: gpu`), or if any explicit `other_parameters.pressureProfile` setting in a stage is truthy on GPU; CUDA-enabled NAMD does not support `pressureProfile`, and the previous silent-strip behavior masked user misconfiguration; the error message names the offending stage's ensemble and points at the two ways to fix the conflict (disable the flag, or use a CPU NAMD build)
- doc: ex16 and ex17 example YAMLs (membrane-embedded HIV-1 gp41 MPER-TM trimer builds) now carry a header comment explaining that their post-equilibration NPAT stage and packaged production run use `pressureProfile`, so they must be run with a CPU NAMD build; the comment names the specific `other_parameters` blocks to strip if a user wants to run them on GPU

## [2.4.12] - 2026-05-13

- bugfix: `PackmolLogParser` failed on logs produced by `packmol` 21.x because the new release reorganized its header in three breaking ways: the `Periodic boundary condition activated:` line became multi-line with indented `Minimum coordinates:` and `Maximum coordinates:` rows, the `PBC Reference box:` line was removed entirely (its information now carried by the min/max rows), and the per-type `Maximum number of GENCAN loops for type:` lines were dropped in favor of the aggregate `for all molecule packing:` value alone; the parser now auto-detects the PBC format by checking whether numeric tokens trail the activation line and reconstructs both `pbc_boxsize` (max-min) and `pbc_reference_box` (concatenated min+max) so downstream consumers see the same metadata schema; logs from packmol 20.15.1 continue to parse correctly via the legacy code path
- bugfix: `PackmolLogParser` initialized each structure's `current_gencan_loops` counter only inside the per-type `Maximum number of GENCAN loops for type:` block, so a `KeyError` was raised on the first GENCAN report for any structure whose type lacked that header line (always the case in packmol 21.x); the counter is now initialized to 0 at structure-creation time
- bugfix: `PackmolLogParser.finalize()` crashed with `'Axes' object is not subscriptable` for single-structure packs because `plt.subplots(1, 1)` returns a bare `Axes` rather than a 1-element array; `squeeze=False` plus `ax.flatten()` makes the indexing uniform

## [2.4.11] - 2026-05-13

- new feature: integrated the VMD `ssrestraints` plugin into the `md` task as a peer of `constraints` and `colvar_specs`; when a top-level `ssrestraints` sub-dict appears under an `md` task, pestifer invokes the plugin against the current state's PSF/PDB to emit a NAMD `extraBonds` file and wires `extraBonds=on` / `extraBondsFile=<basename>.extrabonds` into the generated NAMD config; the schema entry exposes all 13 plugin options (selection text, force constants for protein/NA, `na` mode, `ideal`, H-bond restraint sub-options) with defaults that mirror the plugin's own; the generated extraBonds file participates in the artifact pipeline via a new `NAMDExtraBondsFileArtifact`
- enhancement: `packmol` is now version-validated at config-load time; pestifer probes the resolved binary by sending it an empty stdin and parses the version banner, refusing any version below 20.15.1 (the first release with the `pbc` keyword); when the resolved binary lives inside a conda environment, the error message hints at the AmberTools/`packmol-memgen` shadow problem and recommends installing a standalone packmol outside the conda env; `paths.packmol` may be set to an absolute path in user YAMLs to bypass `PATH` resolution entirely; ex16 and ex17 example YAMLs now demonstrate this pattern
- enhancement: the `PackmolScripter` post-run log scan now flags lines containing `ERROR`, `Could not open file`, or `Could not find` and forces a non-zero return code so callers raise; older packmol could exit 0 while still failing internally, and this surfaces those failures
- bugfix: `[package vcompare [vmdinfo version] X.Y]` returns unreliable results on VMD builds whose version banner has trailing tags or non-numeric components; replaced the comparison in both `pestifer/resources/tcl/macros.tcl` and the macro generator in `pestifer/core/labels.py` with `[lindex [split [vmdinfo version] .] 0] < MAJOR`, an integer-major-version gate that is robust across builds
- doc: explain that the `mdplot` subcommand and task emit per-quantity CSV files alongside plots; explain that `density` is a pestifer-computed timeseries column derived from the cell volume, not a NAMD output; updated the `mdplot` 2.4.x enhancement docs (time-on-x-axis, automatic unit selection)
- test: corrected `_ARCHIVE_PREFIX` in `test_desolvate.py` to match the schema default (`'artifacts'`) so the test reflects what `terminate` actually produces by default; relaxed the per-binary path assertions in `test_config_help_examples` from strict equality to `endswith()` so legitimate per-example `paths.packmol` overrides (such as the absolute paths added to ex16/ex17) do not break the smoke test

## [2.4.10] - 2026-05-11

- enhancement: GPU mode is now auto-detected by comparing `paths.namd3` and `paths.namd3gpu`; if the two paths differ and `namd3gpu` is found in PATH, pestifer switches to GPU mode automatically — no `namd.processor-type` setting required; `processor-type` is retained in the schema for backward compatibility but is ignored; `paths.namd3gpu` defaults to `namd3` (CPU mode when identical); workstation users with separate CPU and GPU-resident NAMD3 builds simply set `paths.namd3gpu` to the GPU binary name/path; the `--gpu` CLI flag now correctly propagates GPU mode to both the config and the live NAMD scripter
- bugfix: GPU-resident NAMD3 config parameter renamed from `GPUResident` to `CUDASOAintegrate` to match current NAMD3 release notes; `outputEnergies: 500` added to the `gpu-resident` schema section (NAMD recommends 500–1000 for GPU-resident mode); `writescript` now suppresses any param key that appears in `gpu-resident` (case-insensitive) so each parameter appears exactly once in the generated config file
- enhancement: multi-GPU GPU-resident launches now include `+pmepes K` to assign fewer PEs to the PME device and reduce load imbalance; K is computed as total PEs minus 8×(N−1) for N GPUs, matching the NAMD3 recommended work distribution examples

## [2.4.9] - 2026-05-11

- enhancement: `mdplot` timeseries plots now use simulation time on the bottom x-axis instead of raw timestep index; units are auto-selected by magnitude (ps when total duration < 1000 ps, ns otherwise); the integer timestep index is shown on the top spine via a secondary x-axis; `NAMDLogParser.finalize()` adds a `dt_fs` constant column to the energy dataframe so that chained runs with different timestep sizes are handled correctly — simulation time is computed in the plot task as `cumsum(dt_fs × ΔTS / 1000)` across the concatenated dataframe, which correctly accumulates time across run segments regardless of whether the timestep changes between segments
- enhancement: `mdplot` timeseries y-axis unit labels now default to the correct physical units for known NAMD output columns (kcal/mol for energy terms, K for temperature, bar for pressure, Å³ for volume) rather than `*`; energy quantities whose maximum absolute value exceeds 1000 kcal/mol are automatically scaled by 1/1000 and labeled `1000 kcal/mol` to avoid unwieldy tick magnitudes

## [2.4.8] - 2026-05-11

- bugfix: `scripts/release.sh` did not update `CITATION.cff`, causing the release workflow to fail when it checked that the citation version matched the release; `release.sh` now bumps both `version` and `date-released` in `CITATION.cff` alongside `pyproject.toml`

## [2.4.7] - 2026-05-11

- bugfix: `pestifer mdplot --timeseries` crashed with `KeyError: 'pressureprofile'` when `--profiles` was not supplied; the `--profiles` argument defaulted to `['pressure']`, causing the task to look for a pressure-profile dataframe that was never loaded; fixed by changing the default to `[]`; also added a guard in the task to filter `profiles` to only those whose dataframe was actually loaded, so the same crash cannot recur even when profiles are explicitly requested but absent from the log files
- bugfix: `pestifer.core.config` imported `Yclept` via `from ycleptic.yclept import Yclept`, a path that no longer exists after ycleptic's 2.0.5 `src/` restructure; updated to `from ycleptic import Yclept` using the new stable public API (requires ycleptic >= 2.0.6)

## [2.4.5] - 2026-05-11

- bugfix: `CharmmParamFile.merge()` incorrectly deduplicated multi-term dihedral parameters; CHARMM sums all dihedral terms for the same atom-type quartet that have different periodicities `n`, so deduplication must key on `(canonical_quartet, n)` rather than the quartet alone; the old code collapsed all terms for a quartet to a single entry, silently discarding all but one periodicity; added five static key methods (`_bond_key`, `_angle_key`, `_dihedral_key`, `_improper_key`, `_nbfix_key`) and rewrote `merge()` to use them; also added 8 unit tests covering bond/dihedral/nonbonded deduplication and multi-term dihedral preservation

## [2.4.4] - 2026-04-28

- bugfix: graft alignment VMD selections produced zero atoms for all grafts in structures where native glycan resids were renumbered by the narrow convention; the N-point alignment code introduced in v2.3.0 used post-renumbering resids from `G.residues` instead of the original shortcode resids (`G.target_root` / `G.target_partners`); since `$m4` is the raw source PDB which has not yet been renumbered at the time graft alignment runs in the Tcl script, the selections must use the original PDB resids; reverted to the pre-v2.3.0 approach of using `G.target_root` and `G.target_partners` directly, generalized to N alignment points
- bugfix: terminus-management patches (`NTER`, `CTER`, `GLYP`, `PROP`, `ACE`, `CNEU`, etc.) and their pestifer undo-counterparts (`XCTR`, `XNTR`, `XGLP`, `XPRP`) recorded in `REMARKS patch` lines of a continuation PSF were being re-propagated into the next psfgen script as generic patches; `CTER`/`NTER` loop-boundary records were emitted as `last`/`first` segment directives (applying to the wrong, segment-terminal residues), while `XCTR`/`XNTR` records were emitted as standalone `patch` commands against residues whose C=O / N-HN bonds had never been removed, producing duplicate bonds and a degenerate `C O C` angle that NAMD cannot parameterize; fixed by filtering all terminus-management and undo patches out of `PSFContents.generic_patches`

## [2.4.3] - 2026-04-28

- bugfix: generic single-residue patches (NTER, CTER, GLUP, ASPP, NNEU, and custom patches) recorded in `REMARKS patch` lines of an input PSF were silently dropped when a `psfgen` task followed a `continuation` task; `PSFContents` now extracts these patches into `generic_patches` and `AsymmetricUnit` propagates them through the objmanager pipeline so they are re-emitted as `first`/`last`/`patch` commands in the new psfgen script
- bugfix: `ring_check` task now operates correctly without an XSC file (vacuum / non-periodic systems); previously the task assumed a periodic box was always available; in non-periodic mode the bounding box is derived from coordinate extents and minimum-image wrapping is skipped
- test: added `tests/unit/test_tasks/test_psfgen.py` — slow regression test that builds a vacuum BPTI PSF/PDB (6pti, mirroring example 1), then loads it via a `continuation` task and re-runs `psfgen` with one non-cysteine mutation, verifying that all three DISU patch REMARK records survive into the stage-2 PSF

## [2.4.2] - 2026-04-28

- bugfix: `PSFContents` used `BaseObjList.get()` to filter atoms by segtype when building the topology-active atom list for `ring_check`; `get()` returns `None` for zero matches and a bare object for exactly one match, causing a `TypeError` on `len()`; replaced with `filter()`, which always returns a `BaseObjList`
- new feature: custom exception hierarchy (`PestiferError`, `PestiferBuildError`) replaces bare `RuntimeError`/`ValueError`/`Exception` raises throughout the task and molecule layers; the CLI top-level handler catches any `PestiferError` and logs a clean error message with `sys.exit(1)` rather than printing a Python traceback; `ring_check` with `delete: none` now raises `PestiferBuildError` instead of silently logging

## [2.4.1] - 2026-04-27

- change: `ring_check` task now terminates the build with an error when any glycan ring is found to be pierced; because glycan residues are covalently bonded to the protein they cannot be deleted, so pestifer logs each offending segment and residue and raises a `PestiferBuildError`; the recommended remedy is to rebuild with a different random seed or add more minimization before the `ring_check` task; non-glycan piercings (e.g. sterol rings) continue to be handled by the `delete` parameter as before

## [2.4.0] - 2026-04-27

- new feature: `setup-vmd` subcommand generates `~/.pestifer/vmd_init.tcl` with the absolute path to pestifer's Tcl root and appends a one-line source hook to `~/.vmdrc`; re-running after upgrading pestifer or switching conda environments keeps VMD's Tcl package paths in sync
- new feature: `align` mod accepts a new `ref_sourceID` field (RCSB PDB ID) as an alternative to `ref_pdb`; the reference structure is downloaded to the working directory automatically, keeping `base_coordinates` and `base_molecule` invariant without requiring a second `fetch` task; downloaded files are registered as pipeline artifacts and swept into the artifacts tarball on cleanup
- new feature: `terminate` task logs a structured system report at the end of every build: output file sizes, topology feature counts (atoms, bonds, angles, dihedrals, impropers, cross-terms from the PSF), and periodic box vectors and lengths from the XSC file
- change: `package.basename` in the `terminate` task now controls only the tarball filename and NAMD configuration script name; coordinate state files (PSF, PDB, COOR, XSC, VEL) retain their terminate-task basename inside the tarball and are no longer renamed
- bugfix: `VMDScripter.writescript()` used substring matching to detect an existing `exit` statement, causing `exit 1` in error-branch Tcl code to suppress injection of the final bare `exit` and leaving VMD hanging after script completion; fixed to require an exact line match
- change: progress bar timer widget now uses `progressbar.utils.len_color` (ANSI-aware) instead of Python's built-in `len` so bars correctly fill to the terminal's right edge regardless of color escape sequences in the timer label

## [2.3.1] - 2026-04-24

- bugfix: `GraftList.assign_residues()` now removes base-structure residues downstream of the outermost receiver (and their links) before applying graft external links; previously, if the base structure carried any sugars beyond the outermost receiver resid, psfgen would apply both the existing glycosidic patch and the new graft external link to the same acceptor oxygen, producing a doubly-bonded OC301 and a `CC3162 OC301 CC3162 CC3161` dihedral that NAMD cannot parameterize

## [2.3.0] - 2026-04-24

- bugfix: six `pdbalias atom` entries for sialic acid (ANE5AC) in `labels.py` incorrectly used `ANE5` as the residue prefix; psfgen applies atom aliases using the post-residue-alias CHARMM name (`ANE5AC`), so these aliases were silently ignored and coordinate assignment for SIA atoms failed; all six entries corrected to `ANE5AC`
- new feature: glycan grafts now support **N-point alignment** with **dihedral pre-conditioning**; the graft shortcode target now accepts any number of alignment residues separated by `#` (e.g., `G_572#573#574`), and the source accepts the same number of index resids (e.g., `4b7i,C_1#2#3-8`); before computing the `measure fit` transformation, pestifer matches each glycosidic dihedral (phi, psi, and omega for 1→6/2→6 linkages) in the source stem to the corresponding angle in the target stub by rotating the outer branch around the appropriate bond axis; this substantially improves donor atom placement when the target already carries a multi-sugar stub

## [2.2.11] - 2026-04-24

- bugfix: `toppar_all36_carb_imlab.str` added to the default CHARMM stream file list; this file contains supplementary carbohydrate parameters (e.g., `CC3162 OC301 CC3162` angle) needed for certain glycan-glycan ether linkages in N-glycan trees; its absence caused NAMD to abort with "UNABLE TO FIND ANGLE PARAMETERS" on glycosylated systems built from multi-donor glycan grafts
- new feature: `ring_check` task now includes `glycan` segments in its default `segtypes` (previously only `lipid`); glycan pyranose rings are now detected and checked for lipid-chain piercing by default; note that for glycan rings pierced by lipid chains, `delete: piercer` is more appropriate than the default `delete: piercee`
- new feature: `mdplot` task now saves PNG output to a subdirectory (default `mdplots/`) specified by the new `output_dir` spec; plots are marked `keep=True` so they survive the artifacts tarball sweep and are immediately visible after the run
- change: `mdplot` default colormap changed from `viridis` to `tab10` for better perceptual separation of discrete traces on multi-trace plots
- change: `mdplot` single-trace plots no longer show a legend even when `legend: True` is set

## [2.2.10] - 2026-04-20

- bugfix: graft source molecules (e.g., a standalone glycan PDB fetched for grafting) are no longer subjected to glycan chainID/segname/resid reassignment; the original chainIDs and resids are preserved so that `activate()` can correctly locate the source segment by its original chainID
- bugfix: graft donor residues now go into a dedicated psfgen segment (e.g., `NG01`) instead of being mixed into the protein segment; the segname is computed from the target chain and existing glycan count, consistent with the base molecule glycan naming convention
- bugfix: graft bond patch segnames now correctly reference the base molecule's psfgen segnames; receiver-side residues are looked up via `res_segname_map`, donor-side residues use the graft's dedicated segname (`graft_segname`), eliminating stale source-PDB chainID references in patch commands
- test: added `test_graft_segname_and_link_segnames` to `tests/unit/test_objs/test_graft.py` covering the full segname/resid pipeline for a glycan graft
- test: added `tests/unit/test_molecule/test_graft_molecule.py` — integration tests using real PDB structures (4b7i chain C glycan grafted onto 1gc1 chain G protein) verifying activation counts, `graft_segname` assignment, and link segname correctness end-to-end through the AsymmetricUnit pipeline

## [2.2.9] - 2026-04-20

- new feature: two glycan resid numbering conventions are now available via `sequence.glycans.numbering`
  - `narrow` (default): glycan blocks start above the maximum protein resid of the parent chain; block size controlled by `sequence.glycans.max_glycan_size` (default 30)
  - `wide`: the root glycan monomer (directly bonded to the protein) gets the protein attachment resid plus `sequence.glycans.wide_shift` (default 3000); subsequent monomers increment by 1 in BFS order from the root (e.g. NAG on ASN 123 → resid 3123, next monomer → 3124)
- bugfix: `--ncpus` CLI option and `namd.ncpus` config key now correctly override NAMD PE count; the processor-info banner reflects the effective count at startup
- bugfix: under SLURM, NAMD (CPU/multicore mode) is now launched as `numactl --localalloc namd3 +p N script.namd`; both `charmrun +p N` and bare `+p N` bypass SLURM's CPU affinity mask, scattering Charm++ worker threads across all NUMA nodes and causing severe startup slowdown due to cross-NUMA memory access during bond-table construction; `numactl --localalloc` preserves Charm++'s thread placement while ensuring each thread's memory allocations land on its local NUMA node
- bugfix: glycan resid block allocation in the `narrow` convention now correctly starts above the parent chain's maximum protein resid rather than from 1, eliminating clashes between glycan and protein resids within the same chainID

## [2.2.7] - 2026-04-19

- new feature: all glycan trees (both polymeric and monomeric) now share the chainID of their attached protein chain, each receives a unique segname (`{chainID}G{k:02d}`), and resids are block-allocated per chain in attachment-resid order so each glycan owns a non-overlapping resid range allowing future monomer additions; controlled by new `sequence.glycans.max_glycan_size` parameter (default 30)
- new feature: link `resid1`/`resid2` fields are synchronized from residue objects after segment building to keep topology consistent through `prune_topology`

## [2.2.6] - 2026-04-19

- bugfix: stale `pstr()` call on `StateInterval` in loop-declashing debug log replaced with valid attribute access
- bugfix: `assign_obj_to_attr` in `BaseObj` now correctly handles multi-match lists — previously used `type(x) != objList` (always `True`) instead of `type(x) != type(objList)`, causing a `ResidueList` to be assigned where a single residue was expected, leading to downstream `AttributeError`
- bugfix: `prune_topology` in `SegmentList` now uses identity (`is`) rather than equality (`==`) when locating and removing residues — after all glycan chainIDs are remapped to their parent protein chain, glycan residues from different source chains (e.g., chain A NAG1 and chain D NAG1) become value-equal but are distinct objects; using `__eq__` caused the wrong segment's residues to be deleted, breaking the glycan patch application
- bugfix: `unregister_chain` is no longer called for pruned glycan segments (e.g., GG03) whose `segname` is not a registered chainID — previously raised `KeyError` when a glycan segment was fully pruned by `prune_topology`
- bugfix: asymmetric-unit builder now detects C-terminal insertions and ensures they are built even when `include_C_termini: False` is set

## [2.2.5] - 2026-04-18

- bugfix: `merge` task now writes injected PSF REMARKS with a leading space, consistent with psfgen's own output format

## [2.2.4] - 2026-04-18

- Exclusion logic expressions now tolerate atom-only attributes (e.g., `altloc`) on objects that lack them, returning `False` rather than raising an error — so `altloc == 'B'` can be used safely in a `psfgen` source `exclude` list

## [2.2.3] - 2026-04-17

- Now generates a custom "minimal" CHARMM parameter file customized for every system
- New `merge` task for merging independent PSF/PDBs into a single system (no coordinate transformations!)

## [2.2.0] - 2026-03-04

- Now uses Feb 26 CHARMMFF release (Jul 24 is still available as an option)

## [2.1.9] - 2026-03-03

- bugfix: placeholder residues (insertions) now inherit correct segnames
- bugfix: psfgen task can now follow a continuation task and add atoms and residues successfully

## [2.1.8] - 2026-02-18

- bugfix: desolvate no longer assumes solvent and non-solvent atoms must have different segnames (but they should)

## [2.1.7] - 2025-11-19

- bugfix: insertion resid incrementation bug fixed

## [2.1.6] - 2025-10-22

- bugfix: glycan declashing

## [2.1.5] - 2025-10-22

- testing `GLYCAN_PENDANT` crotation

## [2.1.4] - 2025-10-12

- multiple bugfixes in `desolvate` and `make-namd-restart` subcommands

## [2.1.3] - 2025-10-10

- multiple bugfixes in `make-pdbrepository`
- dynamic patch for `toppar_all36_lipid_cardiolipin.str` to correct bad IC's
- deleted unnecessary double-counting check for glycosylation-type validation
- added cardiolipin substream to built-in pdbrepository

## [2.1.2] - 2025-09-30

- bugfix for orienting transmembrane proteins before calculating membrane x,y dimensions
- bugfix for grafts onto targets that have H's from sources that do not

## [2.1.1] - 2025-09-23

- bugfixes in `mdplot` (pressure profiles) and `make_membrane_system` (`dum_pdb`)
- updated input configs for examples 16 and 17

## [2.1.0] - 2025-09-19

- bugfixes in cif-to-pdb resid mapping output, production NAMD script writer, and mispelling of rcsb
- new feature: `test_standards` optional subtask added to `terminate` task to enable generation of gold standard results and comparisons to them

## [2.0.3] - 2025-09-15

- `mdplot` now correctly treats `cpu_time` and `wall_time` as running sums over chained MD runs
- bugfixes in `new-system` subcommand and `terminate` task
- bugfix in `bilayer_embed.tcl` for z-translation
- updated `pidibble` dependency to allow for downloading structures from OPM

## [2.0.1] - 2025-09-04
- on-the-fly patching `toppar_all36_prot_modify_res.str` so that `PRES ZNHD` and `PRES ZNHE` are correct
- `validate` task introduced to allow for validation of psfgen-produced PSF and PDB files
- `include` and `exclude` now allow full logical expressions for better control of inclusion and exclusion of atoms and residues from source structures
- `fetch` task now separately responsible for downloading necessary pdb files
- `continuation` task now separately responsible for starting from a given state (PSF/PDB/COOR/XSC)
- `mdplot` now correctly treats `cpu_time` and `wall_time` as running sums over chained MD runs
- complete refactoring to use pydantic's BaseModel

## [1.21.2] - 2025-07-17
- first official Zenodo release [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16051499.svg)](https://doi.org/10.5281/zenodo.16051499)
- fixed bug in user-modifiable segtype handling

## [1.21.1] - 2025-07-15
- Parameterized BMS-806 (resname 83G) from PDB ID 8fad now included

## [1.20.3] - 2025-07-14
- more examples
- upgraded example support

## [1.19.1] - 2025-07-11
- fixed namd restart bugs
- fixed improper sorting of segment indices in bilayer_embed.tcl
- support for adding, renaming, inserting, assigning authorship to, and deleting examples

## [1.18.1] - 2025-07-08
- switched from patches to mutations for protonated histines predicted by `pdb2pqr`
- developer support for adding examples

## [1.17.0] - 2025-07-05
- support for simple nucleic acids added

## [1.16.3] - 2025-07-03
- fixed None-return bug in terminate task
- fixed bug in insulin example to correctly assign NE2-Zn bonds
- fixed bug causing piericing of PRO807 in model-building missing residues in insulin receptor example
- new `pdb2pqr` task to assign protonation states based on `propka3`

## [1.15.1] - 2025-06-22
- fixed the dropped filename extensions in `mdplot`
- updated salt concentration implementation in membrane system builder
- removed spurious dependence on `parmed`
- added pressure profile calculation capabilities to `mdplot` task

## [1.14.1] - 2025-06-18
- Implemented standard CHARMM36 patches as mods
- Added HIV-1 protease 1f7a as a new example (now there are 16!)

## [1.13.3] - 2025-06-17
- Allow config to contain names of otherwise unnamed lipids in CHARMMFF

## [1.13.2] - 2025-06-16
- Example 15 added with ring-check
- changed `make-resi-database` to `make-pdb-collection`
- updated PDB repository to work with a gzipped folder (all lipids)
- converted charmmff contents back into an as-downloaded tarball

## [1.12.1] - 2025-06-04
- added `half_mid_zgap` parameter to ease membrane packing
- added `C6DH` and `C7DH` lipid residue aliases for `C6DHPC` and `C7DHCP`, respectively
- enable inadvertantly unimplemented user selection of number of lipids per leaflet in a patch
- misspelled "paramfiles" in ycleptic basefile (one time)

## [1.12.0] - 2025-06-03
- documentation upgrades
- better progress bar integration with `packmol`
- fixed residue misnumbering for large membranes
- fixed large membrane embedding errors

## [1.11.2] - 2025-05-27
- fixed parsing error converting specstrings to composition dicts
- tests passed after minor updates
- New bilayer build protocol
- New NAMD and packmol on-the-fly log parsing, including in standalone mode with the `follow-namd-log` subcommand

## [1.10.0] - 2025-03-27
- Updated the `pestifer_init` TcL script definition to allow the command `pestifer_init` to replace the clunkier `source [pestifer_init]`
- `cleanup` subcommand added for cleaning up after an aborted run

## [1.9.0] - 2025-03-27
- `transrot` mod enabled for global translations and rotations
- `desolvate` subcommand now can output a single pdb file

## [1.8.3] - 2025-03-19
- fixed bug for growing alpha helices in the N-terminal direction in `crot.tcl`

## [1.8.2] - 2025-03-11
- fixed failure to avoid deprecated namd parameters if using namd2
- fixed hard-coded C-direction brot call to detect whether brot is N- or C-terminal
- added `--gpu` option for the `run` and `run-example` subcommands

## [1.8.1] - 2025-03-10
- removal of tcllib from resources (I never used it anyway)
- bugfix in slow build tests

## [1.8.0] - 2025-03-07
- new source code structure
- mmCIF convention updated
- `mdplot` subcommand added

## [1.7.4] - 2025-02-23
- bugfix: xst trace for `mdplot` fails if empty

## [1.7.2] - 2025-02-12
- bugfix: correct the detection of whether or not a file is a NAMD log file based on the first two Info: records
- bugfix: correct the problem with restarts interfering with mdplot

## [1.7.1] - 2025-02-10
- `make-namd-restart` enhanced with automatic SLURM script updating
- bugfix: GPU-resident NAMD3 cannot do mulitple-gpu runs with constraints
- all lipids in the charmmff topology file `top_all36_lipid.rtf` now have PDBs ready for `packmol`
- `show-resources` subcommand enabled
- `make-namd-restart` subcommand enabled
- `ycleptic` dependency updated to 1.1.0

## [1.6.1] - 2025-06-29
- now able to use namd2 or namd3
- can optionally use GPU-resident namd3

## [1.5.9] - 2025-01-22
- bugfix: log write suppressed inadvertently if progress bars not used

## [1.5.6] - 2025-01-03
- `desolvate` subcommand implemented
- wildcard allowed in pdbalias commands for atom renaming
- temporary fix for dbRes HIS in any mutations to be named HSD

## [1.5.4] - 2024-11-05
- bugfix: incorrect deletion of image seqmods

## [1.5.3] - 2024-09-30
- bugfix: `custom_pdb_path` bug in `bilayer` fixed

## [1.5.2] - 2024-09-24
- glycan graph mistake fixed
- python dependency updated to >=3.12
- `ycleptic` dependency updated to 1.0.7

## [1.4.8] - 2024-09-24
- updated CHARMM lipid PDB files
- updated `ycleptic` to 1.0.6 to enable interactive help and automatic config documentation

## [1.4.7] - 2024-09-18
- `ambertools` dependency removed
- `packmol-memgen` integration removed; now use native `bilayer` task
- `make-resi-database` command added
- CHARMM force field files updated to July 2024
- `salt_con`, `anion`, and `cation` specs for solvate now available
- `pidibble` dependency updated to 1.1.9
- pierced ring detection and remediation via the `ring_check` task
- `restart` task added
- automatic detection of SLURM environment for multi-node MD runs
- `--config-updates` option for `fetch-example` and `run-example` subcommands implemented
- progress bars enabled for NAMD, psfgen, and packmol
- `--kick-ass-banner` option implemented -- check it out!
- `pidibble` dependency updated to 1.1.8
- expanded integration of `packmol-memgen`
- added `fetch-example` subcommand that just copies the respective example YAML file to the CWD
- bugfixes:
- since packmol-memgen sometimes translates the insert, cannot use packmol's input coordinates to psfgen the resulting embedded system

## [1.4.4] - 2024-07-10
- now includes Tcllib 2.0
- bugfixes:
- fixed incorrect charges on the C-terminal CA and HB in the `HEAL` patch

## [1.4.3] - 2024-07-02
- update ambertools version requirement to 23.6; no more packmol-memgen/pdbremix error
- bugfixes:
- change packmol-memgen's weird ion names to be CHARMM-compatible
- allow for N-atom position calculation for residues added to a C-terminus (atom name OT1 vs O)

## [1.4.2] - 2024-06-27
- explicit chain mapping in config file

## [1.4.1] - 2024-05-16
- support for empty TER records

## [1.4.0] - 2024-04-01
- initial `packmol-memgen` integration

## [1.3.9] - 2024-03-04
- added `include_C_termini` boolean to `declash` directives; set to `False` to prevent C-terminal insertions from undergoing automatic declashing

## [1.3.8] - 2024-02-29
- bugfix: spurious code in `pestifer-vmd.tcl`

## [1.3.7] - 2024-02-29
- bugfix: fixed a spurious hard-coded path in `macros.tcl`
- bugfix: `runscript` sources TcL proc files with dependencies in proc files that aren't yet sourced; fixed that
- `alphafold` source directive added to permit download of models from the AlphaFold database by accession code

## [1.3.5] - 2024-02-26
- bugfix: renumbering of author resids in non-protein segments if user adds protein residues by insertion that may conflict
- transferance of atomselect macros from YAML input to any VMD script
- `inittcl` subcommand makes this transfer; needs only to be run one time post-installation

## [1.3.4] - 2024-02-06
- new TcL procs for asymmetric unit generation from non-symmetric assemblies
- `pestifer_init` TcL proc provided in docs for user VMD startup script
- `script` subcommand removed
- syntax of `wheretcl` subcommand expanded

## [1.3.3] - 2024-01-31
- `NAMDLog` class introduced for parsing NAMD2-generated log files
- `mdplot` task for generating plots of various energy-like quantities vs timestep

## [1.3.2] - 2024-01-24
- allow for user-defined links in the config file
- all example builds now have tests in the test suite

## [1.3.1] - 2024-01-12
- bug fixes for cleaving

## [1.3.0] - 2024-01-11
- Support for reading from already-built PSF/PDB systems

## [1.2.9] - 2023-12-19
- improved declashing and domain-swapping

## [1.2.8] - 2023-12-05
- `grafts` for adding glycans
- `cleave` task and `CleavageMod`
- `ModManager` replaces `ModContainer`

## [1.2.5] - 2023-11-28
- `insertion` mod; corrected bug in `brot` tcl procedure

## [1.2.3] - 2023-11-20
- script subcommand handles local scripts
- added `wheretcl` subcommand
- added `script` subcommand (since removed)

## [1.2.0] - 2023-11-16
- split all namd subtasks out; now they are level-1 tasks
- added `manipulate` task

## [1.1.2] - 2023-11-09
- more control over production NAMD2 config generated by the package directive
- position restraints control in minimization and relaxation
- `other_parameters` for any NAMD2 relaxation task

## [1.0.9] - 2023-11-07
- alternate coordinate files and Cfusions
- chain-specific control over building in zero-occupancy residues at N and C termini
- `alpha` crotation for folding a span of residues into an alpha helix

## [1.0.6] - 2023-10-31
- `cif_residue_map_file` generated to report mapping between CIF-residue numbering and author residue numbering
- enhancements to packaging task
- support for topogromacs added

## [1.0.1] - 2023-09-20
- Initial release