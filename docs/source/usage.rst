Usage
=====

Installation of the ``pestifer`` package gives access to the ``pestifer`` command.  The general syntax for invoking ``pestifer`` is

.. code-block:: bash

   $ pestifer <subcommand> <options>

Help with any subcommand can be obtained via

.. code-block:: bash

   $ pestifer <subcommand> --help

General options for all subcommands:

.. code-block:: bash

  --no-banner          turn off the banner
  --loglevel LOGLEVEL  Log level for messages written to diagnostic log (debug|info)
  --diag DIAG          diagnostic log file

There are several subcommands:

``run``: the main subcommand that uses a user's input config file to build a system.

.. code-block:: bash

   $ pestifer run <config.yaml>


``run-example``: there are 18 example systems; to run number four, for example:

.. code-block:: bash
   
   $ pestifer run-example 4

``config-help``: Interactive help in constructing a config file.

.. code-block:: bash

   $ pestifer config-help
   Help on user-provided configuration file format
       Help available for charmmff, psfgen, namd2, title, paths, tasks
   $ pestifer config-help tasks
   Help on user-provided configuration file format
   tasks:
      Specifies the tasks to be performed; each is a dictionary with a
         heading which is a reserved task name
      type: list
      Help available for psfgen, ligate, solvate, relax, terminate
   $ pestifer config-help tasks psfgen
   Help on user-provided configuration file format
   tasks->
   psfgen:
      Parameters controlling initial psfgen run
      type: dict
      Help available for source, mods, minimize, cleanup


Example 1
---------

An example YAML configuration file to build the 6PTI system in the psfgen user manual is below:

.. code-block:: yaml

   title: BPTI
   tasks:
      - psfgen:
            source:
            id: 6pti
      - solvate:
      - relax:
            ensemble: NPT  

This can be run (preferably in a clean directory) via

.. code-block:: bash

   $ pestifer run-example 1

Or, alternatively, pasting that content into a local file ``myconfig.yaml``:

.. code-block:: bash

   $ pestifer run myconfig.yaml

Generally, pestifer is instructed to execute a series of ``tasks``.  The first is a ``psfgen`` task, where the experimental coordinates are used with psfgen to generate an initial topology/coordinate file couple. The ``source`` block declares we will use the 6pti entry of the RCSB, and we will generate the default biological assembly (which is the only one for this file). The second task is solvation and ionization. Finally, the system is relaxed at a temperature of 300 K and a pressure of 1 bar.  By default, there is a minimization in the ``psfgen`` task and the ``solvate`` task.

This run will generate a lot of files, but one that is interesting to look at is ``00-complete.yaml``, which is a "filled-in" version of the config file, showing all defaults:

.. code-block:: yaml

   # Ycleptic v 1.0.3.2 -- Cameron F. Abrams -- cfa22@drexel.edu
   # Dump of complete user config file
   charmmff:
      custom:
         parameters:
           - toppar_water_ions.str
           - toppar_all36_moreions.str
         topologies:
           - toppar_water_ions.str
           - toppar_all36_moreions.str
      standard:
         parameters:
           - par_all36m_prot.prm
           - par_all36_carb.prm
           - par_all36_lipid.prm
           - par_all36_carb.prm
           - par_all36_na.prm
           - par_all36_cgenff.prm
           - stream/carb/toppar_all36_carb_glycopeptide.str
           - stream/prot/toppar_all36_prot_modify_res.str
         topologies:
           - top_all36_prot.rtf
           - top_all35_ethers.rtf
           - top_all36_cgenff.rtf
           - top_all36_lipid.rtf
           - top_all36_carb.rtf
           - top_all36_na.rtf
           - stream/carb/toppar_all36_carb_glycopeptide.str
           - stream/prot/toppar_all36_prot_modify_res.str
   namd2:
      barostat:
         langevinpiston: true
         langevinpistondecay: 100
         langevinpistonperiod: 200
         langevinpistontarget: $pressure
         langevinpistontemp: $temperature
         useflexiblecell: false
         usegrouppressure: true
      generic:
         1-4scaling: 1.0
         cutoff: 10.0
         exclude: scaled1-4
         outputenergies: 100
         pairlistdist: 11.5
         paraTypeCharmm: true
         switchdist: 9.0
         switching: true
         temperature: 300
      solvated:
         PME: true
         fullElectFrequency: 2
         nonbondedFreq: 1
         pmegridspacing: 1.0
         rigidbonds: all
         stepspercycle: 10
         timestep: 2.0
         wrapAll: true
      thermostat:
         langevin: true
         langevinDamping: 5
         langevinHydrogen: false
         langevinTemp: $temperature
      vacuum:
         dielectric: 80
         fullElectFrequency: 2
         nonbondedFreq: 1
         rigidbonds: none
         stepspercycle: 4
         timestep: 1.0
   paths:
      charmrun: /usr/local/bin/charmrun
      namd2: /usr/local/bin/namd2
      vmd: /usr/local/bin/vmd
   psfgen:
      aliases:
        - atom ILE CD1 CD
        - atom BGLCNA C7 C
        - atom BGLCNA O7 O
        - atom BGLCNA C8 CT
        - atom BGLCNA N2 N
        - atom ANE5 C10 C
        - atom ANE5 C11 CT
        - atom ANE5 N5 N
        - atom ANE5 O1A O11
        - atom ANE5 O1B O12
        - atom ANE5 O10 O
        - atom VCG C01 C1
        - atom VCG C01 C1
        - atom VCG C02 C2
        - atom VCG C03 C3
        - atom VCG C04 C4
        - atom VCG C05 C5
        - atom VCG C06 C6
        - atom VCG C07 C7
        - atom VCG C08 C8
        - atom VCG C09 C9
        - atom TIP3 O OH2
        - residue HIS HSD
        - residue PO4 H2PO4
        - residue MAN AMAN
        - residue BMA BMAN
        - residue NAG BGLCNA
        - residue FUC AFUC
        - residue GAL BGAL
        - residue ANE5 ANE5AC
        - residue SIA ANE5AC
        - residue EIC LIN
        - residue HOH TIP3
        - residue ZN ZN2
        - residue CL CLA
      segtypes:
         glycan:
            resnames:
              - BMA
              - FUC
              - GAL
              - MAN
              - NAG
              - SIA
              - ANE5
         ion:
            resnames:
              - LIT
              - SOD
              - MG
              - POT
              - CAL
              - RUB
              - CES
              - BAR
              - ZN
              - CAD
              - CL
              - SO4
              - PO4
         ligand:
            resnames:
              - EIC
              - VCG
              - 83G
         other:
            resnames: []
         protein:
            invrescodes:
               A: ALA
               C: CYS
               D: ASP
               E: GLU
               F: PHE
               G: GLY
               H: HSE
               I: ILE
               K: LYS
               L: LEU
               M: MET
               N: ASN
               P: PRO
               Q: GLN
               R: ARG
               S: SER
               T: THR
               V: VAL
               W: TRP
               Y: TYR
            rescodes:
               ALA: A
               ARG: R
               ASN: N
               ASP: D
               CYS: C
               GLN: Q
               GLU: E
               GLY: G
               HSE: H
               ILE: I
               LEU: L
               LYS: K
               MET: M
               PHE: F
               PRO: P
               SER: S
               THR: T
               TRP: W
               TYR: Y
               VAL: V
            resnames:
              - ALA
              - ARG
              - ASN
              - ASP
              - CYS
              - GLN
              - GLU
              - GLY
              - HIS
              - HSD
              - HSE
              - ILE
              - LEU
              - LYS
              - MET
              - PHE
              - PRO
              - SER
              - THR
              - TRP
              - TYR
              - VAL
         water:
            resnames:
              - HOH
   tasks:
      - psfgen:
         cleanup: true
         minimize:
            dcdfreq: 100
            nminsteps: 1000
         mods:
            deletions: []
            mutations: []
            ssbonds: []
            ssbondsdelete: []
            substitutions: []
         source:
            biological_assembly: 0
            exclude: {}
            file_format: PDB
            id: 6pti
            sequence:
            fix_conflicts: true
            fix_engineered_mutations: true
            include_terminal_loops: false
            loops:
               declash:
                  maxcycles: 20
               min_loop_length: 4
               sac_res_name: GLY
      - solvate:
         minimize:
            dcdfreq: 100
            nminsteps: 1000
         pad: 10
      - relax:
         dcdfreq: 100
         ensemble: NPT
         nminsteps: 0
         nsteps: 1000
         pressure: 1
         temperature: 300
         xstfreq: 100
   title: BPTI

That's a lot of stuff!  But it shows you everything that Pestifer needs and does to use psfgen and namd2 to generate this system.  Further updates to the documentation will explain in detail.  For now, I recommend running a few of the 18 examples to showcase some of Pestifer's capabilities.  `Ycleptic <https://pypi.org/project/ycleptic/>`_ is a package I developed for generalizing the use of YAML-format configuration files; using a "base" configuration owned by Pestifer, Ycleptic knows how to interpret a user config file to assign defaults, report errors, etc.

The import output of this build are the PSF/PDB/COOR/VEL/XSC files needed to (re)start namd2; by default, these are ``my_system.pdb`` etc.

.. code-block:: bash

   $ ls my_system*
   my_system.coor  my_system.pdb  my_system.psf  my_system.vel  my_system.xsc