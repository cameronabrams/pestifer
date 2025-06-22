# Pestifer
> NAMD Input System Builder

Pestifer is a fully automated simulation-ready MD system builder, requiring as inputs only biomolecular structures (e.g., PDB IDs, PDB files, mmCIF files, alphafold IDs) and a handful of customization parameters, to generate NAMD-compatible input files (PSF, PDB, and xsc).  It is basically a highly functionalized front end for VMD's `psfgen` utility.  It also has a few handy subcommands for working with NAMD output.

## Installation

```bash
pip install pestifer
```

Once installed, the user has access to the main `pestifer` command. 

Pestifer also requires access to the following executables:

1. `namd2` and `charmrun`
2. `vmd` and `catdcd`
3. `packmol`

By default, pestifer looks for these commands in your path.  Specific paths for these can be stipulated in the `paths` directive of your input file.

Pestifer **includes** the [July 2024 Charmm36 force field](https://mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/toppar_c36_jul24.tgz).

## Documentation

Please visit [readthedocs](https://pestifer.readthedocs.io/en/latest) for full documentation.

## Release History
* 1.15.0
    * updated salt concentration implementation in membrane system builder
    * removed spurious dependence on `parmed`
    * added pressure profile calculation capabilities to `mdplot` task
* 1.14.1
    * new docs
* 1.14.0
    * Implemented standard CHARMM36 patches as mods
    * Added HIV-1 protease 1f7a as a new example (now there are 16!)
* 1.13.3
    * Allow config to contain names of otherwise unnamed lipids in CHARMMFF
* 1.13.2
    * Example 15 added with ring-check
* 1.13.1
    * changed `make-resi-database` to `make-pdb-collection`
    * updated PDB repository to work with a gzipped folder (all lipids)
* 1.13.0
    * converted charmmff contents back into an as-downloaded tarball
* 1.12.1
    * added `half_mid_zgap` parameter to ease membrane packing
    * added `C6DH` and `C7DH` lipid residue aliases for `C6DHPC` and `C7DHCP`, respectively
    * enable inadvertantly unimplemented user selection of number of lipids per leaflet in a patch
    * misspelled "paramfiles" in ycleptic basefile (one time)
* 1.12.0
    * documentation upgrades
    * better progress bar integration with `packmol`
* 1.11.3
    * fixed residue misnumbering for large membranes
    * fixed large membrane embedding errors
* 1.11.2
    * fixed parsing error converting specstrings to composition dicts
* 1.11.1
    * tests passed after minor updates
* 1.11.0
    * New bilayer build protocol
    * New NAMD and packmol on-the-fly log parsing, including in standalone mode with the `follow-namd-log` subcommand
* 1.10.0
    * Updated the `pestifer_init` TcL script definition to allow the command `pestifer_init` to replace the clunkier `source [pestifer_init]`
    * `cleanup` subcommand added for cleaning up after an aborted run
* 1.9.0
    * `transrot` mod enabled for global translations and rotations
    * `desolvate` subcommand now can output a single pdb file
* 1.8.3
    * fixed bug for growing alpha helices in the N-terminal direction in `crot.tcl`
* 1.8.2
    * fixed failure to avoid deprecated namd parameters if using namd2
    * fixed hard-coded C-direction brot call to detect whether brot is N- or C-terminal
    * added `--gpu` option for the `run` and `run-example` subcommands
* 1.8.1
    * removal of tcllib from resources (I never used it anyway)
    * bugfix in slow build tests
* 1.8.0
    * new source code structure
    * mmCIF convention updated
    * `mdplot` subcommand added
* 1.7.4
    * bugfix:
        * xst trace for `mdplot` fails if empty
* 1.7.2
    * bugfix: 
        * correct the detection of whether or not a file is a NAMD log file based on the first two Info: records
        * correct the problem with restarts interfering with mdplot
* 1.7.1
    * `make-namd-restart` enhanced with automatic SLURM script updating
* 1.6.4a2
    * bugfix: GPU-resident NAMD3 cannot do mulitple-gpu runs with constraints
* 1.6.4a0
    * all lipids in the charmmff topology file `top_all36_lipid.rtf` now have PDBs ready for `packmol`
    * `show-resources` subcommand enabled
    * `make-namd-restart` subcommand enabled
    * `ycleptic` dependency updated to 1.1.0
* 1.6.1
    * now able to use namd2 or namd3
    * can optionally use GPU-resident namd3
* 1.5.9
    * bugfix: log write suppressed inadvertently if progress bars not used
* 1.5.6
    * `desolvate` subcommand implemented
* 1.5.5
    * wildcard allowed in pdbalias commands for atom renaming
    * temporary fix for dbRes HIS in any mutations to be named HSD
* 1.5.4
    * bugfix: incorrect deletion of image seqmods
* 1.5.3
    * bugfix: `custom_pdb_path` bug in `bilayer` fixed
* 1.5.2
    * glycan graph mistake fixed
    * python dependency updated to >=3.12
    * `ycleptic` dependency updated to 1.0.7
* 1.4.8
    * updated CHARMM lipid PDB files
    * updated `ycleptic` to 1.0.6 to enable interactive help and automatic config documentation
* 1.4.7
    * `ambertools` dependency removed
    * `packmol-memgen` integration removed; now use native `bilayer` task
    * `make-resi-database` command added
    * CHARMM force field files updated to July 2024
    * `salt_con`, `anion`, and `cation` specs for solvate now available
    * `pidibble` dependency updated to 1.1.9
* 1.4.6:
    * pierced ring detection and remediation via the `ring_check` task
    * `restart` task added
    * automatic detection of SLURM environment for multi-node MD runs
    * `--config-updates` option for `fetch-example` and `run-example` subcommands implemented
    * progress bars enabled for NAMD, psfgen, and packmol
    * `--kick-ass-banner` option implemented -- check it out!
    * `pidibble` dependency updated to 1.1.8
    * expanded integration of `packmol-memgen`
* 1.4.5:
    * added `fetch-example` subcommand that just copies the respective example YAML file to the CWD
    * bugfixes:
      * since packmol-memgen sometimes translates the insert, cannot use packmol's input coordinates to psfgen the resulting embedded system
* 1.4.4:
    * now includes Tcllib 2.0
    * bugfixes:
      * fixed incorrect charges on the C-terminal CA and HB in the `HEAL` patch
* 1.4.3:
    * update ambertools version requirement to 23.6; no more packmol-memgen/pdbremix error
    * bugfixes: 
      * change packmol-memgen's weird ion names to be CHARMM-compatible
      * allow for N-atom position calculation for residues added to a C-terminus (atom name OT1 vs O)
* 1.4.2:
    * explicit chain mapping in config file
* 1.4.1:
    * support for empty TER records
* 1.4.0:
    * initial `packmol-memgen` integration
* 1.3.9:
    * added `include_C_termini` boolean to `declash` directives; set to `False` to prevent C-terminal insertions from undergoing automatic declashing
* 1.3.8
    * bugfix: spurious code in `pestifer-vmd.tcl`
* 1.3.7
    * bugfix: fixed a spurious hard-coded path in `macros.tcl`
* 1.3.6
    * bugfix:
      * `runscript` sources TcL proc files with dependencies in proc files that aren't yet sourced; fixed that
    * `alphafold` source directive added to permit download of models from the AlphaFold database by accession code
* 1.3.5
    * bugfix:
      * renumbering of author resids in non-protein segments if user adds protein residues by insertion that may conflict
    * transferance of atomselect macros from YAML input to any VMD script
    * `inittcl` subcommand makes this transfer; needs only to be run one time post-installation
* 1.3.4
    * new TcL procs for asymmetric unit generation from non-symmetric assemblies
    * `pestifer_init` TcL proc provided in docs for user VMD startup script
    * `script` subcommand removed
    * syntax of `wheretcl` subcommand expanded
* 1.3.3
    * `NAMDLog` class introduced for parsing NAMD2-generated log files
    * `mdplot` task for generating plots of various energy-like quantities vs timestep
* 1.3.2
    * allow for user-defined links in the config file
    * all example builds now have tests in the test suite
* 1.3.1
    * bug fixes for cleaving
* 1.3.0
    * Support for reading from already-built PSF/PDB systems
* 1.2.9
    * improved declashing and domain-swapping
* 1.2.8
    * `grafts` for adding glycans
* 1.2.7
    * `cleave` task and `CleavageMod`
* 1.2.6
    * `ModManager` replaces `ModContainer`
* 1.2.5
    * `insertion` mod; corrected bug in `brot` tcl procedure
* 1.2.3
    * script subcommand handles local scripts
* 1.2.2
    * added `wheretcl` subcommand
* 1.2.1
    * added `script` subcommand
* 1.2.0
    * split all namd subtasks out; now they are level-1 tasks
* 1.1.3
    * added `manipulate` task
* 1.1.2
    * more control over production NAMD2 config generated by the package directive
* 1.1.1
    * position restraints control in minimization and relaxation
* 1.1.0
    * `other_parameters` for any NAMD2 relaxation task
* 1.0.9
    * alternate coordinate files and Cfusions
* 1.0.8
    * chain-specific control over building in zero-occupancy residues at N and C termini
* 1.0.7
    * `alpha` crotation for folding a span of residues into an alpha helix
* 1.0.6
    * `cif_residue_map_file` generated to report mapping between CIF-residue numbering and author residue numbering
* 1.0.5
    * enhancements to packaging task
* 1.0.4
    * support for topogromacs added
* 1.0.1
    * Initial release

## Meta

Cameron F. Abrams

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

## Contributing

1. Fork it (<https://github.com/cameronabrams/pestifer/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request

