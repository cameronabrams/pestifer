# Changelog

Pestifer follows [Semantic Versioning](https://semver.org/) and documents changes below.

## [Unreleased]

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