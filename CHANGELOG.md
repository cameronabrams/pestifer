# Changelog

Pestifer follows [Semantic Versioning](https://semver.org/) and documents changes below.

## [UNRELEASED]
### Added
### Changed

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