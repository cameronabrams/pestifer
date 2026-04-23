Changelog
=========

Pestifer follows `Semantic Versioning <https://semver.org/>`__ and
documents changes below.

[2.2.10] - 2026-04-23
---------------------

-  New ``align`` coormod for the ``manipulate`` task: rigidly aligns the
   pipeline system to a reference PDB (or PSF+PDB) by computing the
   least-RMSD rotation and translation via VMD's ``measure fit``.
   Supports independent ``mobile_sel``, ``ref_sel``, and ``apply_to``
   VMD atomselect strings; a congruency check guards against mismatched
   atom counts before the fit is attempted.
-  New ``transfer_coords`` coormod for the ``manipulate`` task: copies
   atomic coordinates from a selected subset of a donor PDB directly onto
   matching atoms of the pipeline system.  An optional pre-alignment step
   (``align_donor_sel`` + ``align_mobile_sel``) rigidly fits the entire
   donor to the pipeline before the transfer, accommodating donors prepared
   in a different coordinate frame.  Both alignment and transfer selections
   are validated for atom-count congruency.
-  ``merge`` task now auto-detects and prepends the pipeline state when
   it is available, eliminating the need for an explicit flag.  An
   optional ``pipeline_system`` key accepts ``segname_map`` and
   ``chainID_map`` sub-keys for renaming the pipelined system's segments
   and chains before the merge.
-  ``merge`` task now supports ``chainID_map`` per input system to
   relabel PDB chain IDs before merging.

.. _section-0:

[2.2.4] - 2026-04-18
--------------------

-  Exclusion logic expressions now tolerate atom-only attributes (e.g.,
   ``altloc``) on objects that lack them, returning ``False`` rather
   than raising an error — so ``altloc == 'B'`` can be used safely in a
   ``psfgen`` source ``exclude`` list

.. _section-1:

[2.2.3] - 2026-04-17
--------------------

-  Now generates a custom “minimal” CHARMM parameter file customized for
   every system
-  New ``merge`` task for merging independent PSF/PDBs into a single
   system (no coordinate transformations!)

.. _section-2:

[2.2.0] - 2026-03-04
--------------------

-  Now uses Feb 26 CHARMMFF release (Jul 24 is still available as an
   option)

.. _section-3:

[2.1.9] - 2026-03-03
--------------------

-  bugfix: placeholder residues (insertions) now inherit correct
   segnames
-  bugfix: psfgen task can now follow a continuation task and add atoms
   and residues successfully

.. _section-4:

[2.1.8] - 2026-02-18
--------------------

-  bugfix: desolvate no longer assumes solvent and non-solvent atoms
   must have different segnames (but they should)

.. _section-5:

[2.1.7] - 2025-11-19
--------------------

-  bugfix: insertion resid incrementation bug fixed

.. _section-6:

[2.1.6] - 2025-10-22
--------------------

-  bugfix: glycan declashing

.. _section-7:

[2.1.5] - 2025-10-22
--------------------

-  testing ``GLYCAN_PENDANT`` crotation

.. _section-8:

[2.1.4] - 2025-10-12
--------------------

-  multiple bugfixes in ``desolvate`` and ``make-namd-restart``
   subcommands

.. _section-9:

[2.1.3] - 2025-10-10
--------------------

-  multiple bugfixes in ``make-pdbrepository``
-  dynamic patch for ``toppar_all36_lipid_cardiolipin.str`` to correct
   bad IC’s
-  deleted unnecessary double-counting check for glycosylation-type
   validation
-  added cardiolipin substream to built-in pdbrepository

.. _section-10:

[2.1.2] - 2025-09-30
--------------------

-  bugfix for orienting transmembrane proteins before calculating
   membrane x,y dimensions
-  bugfix for grafts onto targets that have H’s from sources that do not

.. _section-11:

[2.1.1] - 2025-09-23
--------------------

-  bugfixes in ``mdplot`` (pressure profiles) and
   ``make_membrane_system`` (``dum_pdb``)
-  updated input configs for examples 16 and 17

.. _section-12:

[2.1.0] - 2025-09-19
--------------------

-  bugfixes in cif-to-pdb resid mapping output, production NAMD script
   writer, and mispelling of rcsb
-  new feature: ``test_standards`` optional subtask added to
   ``terminate`` task to enable generation of gold standard results and
   comparisons to them

.. _section-13:

[2.0.3] - 2025-09-15
--------------------

-  ``mdplot`` now correctly treats ``cpu_time`` and ``wall_time`` as
   running sums over chained MD runs
-  bugfixes in ``new-system`` subcommand and ``terminate`` task
-  bugfix in ``bilayer_embed.tcl`` for z-translation
-  updated ``pidibble`` dependency to allow for downloading structures
   from OPM

.. _section-14:

[2.0.1] - 2025-09-04
--------------------

-  on-the-fly patching ``toppar_all36_prot_modify_res.str`` so that
   ``PRES ZNHD`` and ``PRES ZNHE`` are correct
-  ``validate`` task introduced to allow for validation of
   psfgen-produced PSF and PDB files
-  ``include`` and ``exclude`` now allow full logical expressions for
   better control of inclusion and exclusion of atoms and residues from
   source structures
-  ``fetch`` task now separately responsible for downloading necessary
   pdb files
-  ``continuation`` task now separately responsible for starting from a
   given state (PSF/PDB/COOR/XSC)
-  ``mdplot`` now correctly treats ``cpu_time`` and ``wall_time`` as
   running sums over chained MD runs
-  complete refactoring to use pydantic’s BaseModel

.. _section-15:

[1.21.2] - 2025-07-17
---------------------

-  first official Zenodo release |DOI|
-  fixed bug in user-modifiable segtype handling

.. _section-16:

[1.21.1] - 2025-07-15
---------------------

-  Parameterized BMS-806 (resname 83G) from PDB ID 8fad now included

.. _section-17:

[1.20.3] - 2025-07-14
---------------------

-  more examples
-  upgraded example support

.. _section-18:

[1.19.1] - 2025-07-11
---------------------

-  fixed namd restart bugs
-  fixed improper sorting of segment indices in bilayer_embed.tcl
-  support for adding, renaming, inserting, assigning authorship to, and
   deleting examples

.. _section-19:

[1.18.1] - 2025-07-08
---------------------

-  switched from patches to mutations for protonated histines predicted
   by ``pdb2pqr``
-  developer support for adding examples

.. _section-20:

[1.17.0] - 2025-07-05
---------------------

-  support for simple nucleic acids added

.. _section-21:

[1.16.3] - 2025-07-03
---------------------

-  fixed None-return bug in terminate task
-  fixed bug in insulin example to correctly assign NE2-Zn bonds
-  fixed bug causing piericing of PRO807 in model-building missing
   residues in insulin receptor example
-  new ``pdb2pqr`` task to assign protonation states based on
   ``propka3``

.. _section-22:

[1.15.1] - 2025-06-22
---------------------

-  fixed the dropped filename extensions in ``mdplot``
-  updated salt concentration implementation in membrane system builder
-  removed spurious dependence on ``parmed``
-  added pressure profile calculation capabilities to ``mdplot`` task

.. _section-23:

[1.14.1] - 2025-06-18
---------------------

-  Implemented standard CHARMM36 patches as mods
-  Added HIV-1 protease 1f7a as a new example (now there are 16!)

.. _section-24:

[1.13.3] - 2025-06-17
---------------------

-  Allow config to contain names of otherwise unnamed lipids in CHARMMFF

.. _section-25:

[1.13.2] - 2025-06-16
---------------------

-  Example 15 added with ring-check
-  changed ``make-resi-database`` to ``make-pdb-collection``
-  updated PDB repository to work with a gzipped folder (all lipids)
-  converted charmmff contents back into an as-downloaded tarball

.. _section-26:

[1.12.1] - 2025-06-04
---------------------

-  added ``half_mid_zgap`` parameter to ease membrane packing
-  added ``C6DH`` and ``C7DH`` lipid residue aliases for ``C6DHPC`` and
   ``C7DHCP``, respectively
-  enable inadvertantly unimplemented user selection of number of lipids
   per leaflet in a patch
-  misspelled “paramfiles” in ycleptic basefile (one time)

.. _section-27:

[1.12.0] - 2025-06-03
---------------------

-  documentation upgrades
-  better progress bar integration with ``packmol``
-  fixed residue misnumbering for large membranes
-  fixed large membrane embedding errors

.. _section-28:

[1.11.2] - 2025-05-27
---------------------

-  fixed parsing error converting specstrings to composition dicts
-  tests passed after minor updates
-  New bilayer build protocol
-  New NAMD and packmol on-the-fly log parsing, including in standalone
   mode with the ``follow-namd-log`` subcommand

.. _section-29:

[1.10.0] - 2025-03-27
---------------------

-  Updated the ``pestifer_init`` TcL script definition to allow the
   command ``pestifer_init`` to replace the clunkier
   ``source [pestifer_init]``
-  ``cleanup`` subcommand added for cleaning up after an aborted run

.. _section-30:

[1.9.0] - 2025-03-27
--------------------

-  ``transrot`` mod enabled for global translations and rotations
-  ``desolvate`` subcommand now can output a single pdb file

.. _section-31:

[1.8.3] - 2025-03-19
--------------------

-  fixed bug for growing alpha helices in the N-terminal direction in
   ``crot.tcl``

.. _section-32:

[1.8.2] - 2025-03-11
--------------------

-  fixed failure to avoid deprecated namd parameters if using namd2
-  fixed hard-coded C-direction brot call to detect whether brot is N-
   or C-terminal
-  added ``--gpu`` option for the ``run`` and ``run-example``
   subcommands

.. _section-33:

[1.8.1] - 2025-03-10
--------------------

-  removal of tcllib from resources (I never used it anyway)
-  bugfix in slow build tests

.. _section-34:

[1.8.0] - 2025-03-07
--------------------

-  new source code structure
-  mmCIF convention updated
-  ``mdplot`` subcommand added

.. _section-35:

[1.7.4] - 2025-02-23
--------------------

-  bugfix: xst trace for ``mdplot`` fails if empty

.. _section-36:

[1.7.2] - 2025-02-12
--------------------

-  bugfix: correct the detection of whether or not a file is a NAMD log
   file based on the first two Info: records
-  bugfix: correct the problem with restarts interfering with mdplot

.. _section-37:

[1.7.1] - 2025-02-10
--------------------

-  ``make-namd-restart`` enhanced with automatic SLURM script updating
-  bugfix: GPU-resident NAMD3 cannot do mulitple-gpu runs with
   constraints
-  all lipids in the charmmff topology file ``top_all36_lipid.rtf`` now
   have PDBs ready for ``packmol``
-  ``show-resources`` subcommand enabled
-  ``make-namd-restart`` subcommand enabled
-  ``ycleptic`` dependency updated to 1.1.0

.. _section-38:

[1.6.1] - 2025-06-29
--------------------

-  now able to use namd2 or namd3
-  can optionally use GPU-resident namd3

.. _section-39:

[1.5.9] - 2025-01-22
--------------------

-  bugfix: log write suppressed inadvertently if progress bars not used

.. _section-40:

[1.5.6] - 2025-01-03
--------------------

-  ``desolvate`` subcommand implemented
-  wildcard allowed in pdbalias commands for atom renaming
-  temporary fix for dbRes HIS in any mutations to be named HSD

.. _section-41:

[1.5.4] - 2024-11-05
--------------------

-  bugfix: incorrect deletion of image seqmods

.. _section-42:

[1.5.3] - 2024-09-30
--------------------

-  bugfix: ``custom_pdb_path`` bug in ``bilayer`` fixed

.. _section-43:

[1.5.2] - 2024-09-24
--------------------

-  glycan graph mistake fixed
-  python dependency updated to >=3.12
-  ``ycleptic`` dependency updated to 1.0.7

.. _section-44:

[1.4.8] - 2024-09-24
--------------------

-  updated CHARMM lipid PDB files
-  updated ``ycleptic`` to 1.0.6 to enable interactive help and
   automatic config documentation

.. _section-45:

[1.4.7] - 2024-09-18
--------------------

-  ``ambertools`` dependency removed
-  ``packmol-memgen`` integration removed; now use native ``bilayer``
   task
-  ``make-resi-database`` command added
-  CHARMM force field files updated to July 2024
-  ``salt_con``, ``anion``, and ``cation`` specs for solvate now
   available
-  ``pidibble`` dependency updated to 1.1.9
-  pierced ring detection and remediation via the ``ring_check`` task
-  ``restart`` task added
-  automatic detection of SLURM environment for multi-node MD runs
-  ``--config-updates`` option for ``fetch-example`` and ``run-example``
   subcommands implemented
-  progress bars enabled for NAMD, psfgen, and packmol
-  ``--kick-ass-banner`` option implemented – check it out!
-  ``pidibble`` dependency updated to 1.1.8
-  expanded integration of ``packmol-memgen``
-  added ``fetch-example`` subcommand that just copies the respective
   example YAML file to the CWD
-  bugfixes:
-  since packmol-memgen sometimes translates the insert, cannot use
   packmol’s input coordinates to psfgen the resulting embedded system

.. _section-46:

[1.4.4] - 2024-07-10
--------------------

-  now includes Tcllib 2.0
-  bugfixes:
-  fixed incorrect charges on the C-terminal CA and HB in the ``HEAL``
   patch

.. _section-47:

[1.4.3] - 2024-07-02
--------------------

-  update ambertools version requirement to 23.6; no more
   packmol-memgen/pdbremix error
-  bugfixes:
-  change packmol-memgen’s weird ion names to be CHARMM-compatible
-  allow for N-atom position calculation for residues added to a
   C-terminus (atom name OT1 vs O)

.. _section-48:

[1.4.2] - 2024-06-27
--------------------

-  explicit chain mapping in config file

.. _section-49:

[1.4.1] - 2024-05-16
--------------------

-  support for empty TER records

.. _section-50:

[1.4.0] - 2024-04-01
--------------------

-  initial ``packmol-memgen`` integration

.. _section-51:

[1.3.9] - 2024-03-04
--------------------

-  added ``include_C_termini`` boolean to ``declash`` directives; set to
   ``False`` to prevent C-terminal insertions from undergoing automatic
   declashing

.. _section-52:

[1.3.8] - 2024-02-29
--------------------

-  bugfix: spurious code in ``pestifer-vmd.tcl``

.. _section-53:

[1.3.7] - 2024-02-29
--------------------

-  bugfix: fixed a spurious hard-coded path in ``macros.tcl``
-  bugfix: ``runscript`` sources TcL proc files with dependencies in
   proc files that aren’t yet sourced; fixed that
-  ``alphafold`` source directive added to permit download of models
   from the AlphaFold database by accession code

.. _section-54:

[1.3.5] - 2024-02-26
--------------------

-  bugfix: renumbering of author resids in non-protein segments if user
   adds protein residues by insertion that may conflict
-  transferance of atomselect macros from YAML input to any VMD script
-  ``inittcl`` subcommand makes this transfer; needs only to be run one
   time post-installation

.. _section-55:

[1.3.4] - 2024-02-06
--------------------

-  new TcL procs for asymmetric unit generation from non-symmetric
   assemblies
-  ``pestifer_init`` TcL proc provided in docs for user VMD startup
   script
-  ``script`` subcommand removed
-  syntax of ``wheretcl`` subcommand expanded

.. _section-56:

[1.3.3] - 2024-01-31
--------------------

-  ``NAMDLog`` class introduced for parsing NAMD2-generated log files
-  ``mdplot`` task for generating plots of various energy-like
   quantities vs timestep

.. _section-57:

[1.3.2] - 2024-01-24
--------------------

-  allow for user-defined links in the config file
-  all example builds now have tests in the test suite

.. _section-58:

[1.3.1] - 2024-01-12
--------------------

-  bug fixes for cleaving

.. _section-59:

[1.3.0] - 2024-01-11
--------------------

-  Support for reading from already-built PSF/PDB systems

.. _section-60:

[1.2.9] - 2023-12-19
--------------------

-  improved declashing and domain-swapping

.. _section-61:

[1.2.8] - 2023-12-05
--------------------

-  ``grafts`` for adding glycans
-  ``cleave`` task and ``CleavageMod``
-  ``ModManager`` replaces ``ModContainer``

.. _section-62:

[1.2.5] - 2023-11-28
--------------------

-  ``insertion`` mod; corrected bug in ``brot`` tcl procedure

.. _section-63:

[1.2.3] - 2023-11-20
--------------------

-  script subcommand handles local scripts
-  added ``wheretcl`` subcommand
-  added ``script`` subcommand (since removed)

.. _section-64:

[1.2.0] - 2023-11-16
--------------------

-  split all namd subtasks out; now they are level-1 tasks
-  added ``manipulate`` task

.. _section-65:

[1.1.2] - 2023-11-09
--------------------

-  more control over production NAMD2 config generated by the package
   directive
-  position restraints control in minimization and relaxation
-  ``other_parameters`` for any NAMD2 relaxation task

.. _section-66:

[1.0.9] - 2023-11-07
--------------------

-  alternate coordinate files and Cfusions
-  chain-specific control over building in zero-occupancy residues at N
   and C termini
-  ``alpha`` crotation for folding a span of residues into an alpha
   helix

.. _section-67:

[1.0.6] - 2023-10-31
--------------------

-  ``cif_residue_map_file`` generated to report mapping between
   CIF-residue numbering and author residue numbering
-  enhancements to packaging task
-  support for topogromacs added

.. _section-68:

[1.0.1] - 2023-09-20
--------------------

-  Initial release

.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.16051499.svg
   :target: https://doi.org/10.5281/zenodo.16051499
