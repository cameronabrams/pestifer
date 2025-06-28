make-pdb-collection
-------------------

A PDB collection is a set of representative PDB files for small molecules, such as lipids.  Collections are associated with CHARMMFF "streams", and one more more collections comprise a PDB "repository".  A PDB collection must be a directory whose name is its ID, and whose contents are either standalone PDB files with the naming conventions ``<RESI>.pdb``, or subdirectories named after the RESI, each of which contains a set of PDB files that represent different configurations of the molecule. Pestifer's built-in PDB repository was constructed from selected residues in the lipid, water, and ion streams of the CHARMM36 force field.

The built-in PDB collection
+++++++++++++++++++++++++++

You can see what RESIs are included in pestifer using the ``show-resources`` subcommand:

.. code-block:: bash

    $ pestifer show-resources --charmmff pdb
    ---------------------------------------------------------------------------
    PDB Collections:
    PDBCollection(registered_at=2, streamID=water_ions, path=water_ions.tgz, contains 12 resnames)
    BAR,    CAL,    CD2,    CES,    CLA,    LIT,     MG, 
    POT,    RUB,    SOD,   TIP3,    ZN2
    PDBCollection(registered_at=1, streamID=lipid, path=lipid.tgz, contains 130 resnames)
    23SM,    ASM,    BSM, C6DHPC, CER160, CER180, CER181, 
    CER2, CER200, CER220, CER240, CER241,  CER3E,   CHL1, 
    CHM1,   CHSD,   CHSP,   DAPA,   DAPC,   DAPE,   DAPG, 
    DAPS,   DCPC,  DDOPC,  DDOPE,  DDOPS,   DDPC,   DEPA, 
    DEPC,   DEPE,   DEPG,   DEPS,   DGPA,   DGPC,   DGPE, 
    DGPG,   DGPS,   DIPA,   DLPA,   DLPC,   DLPE,   DLPG, 
    DLPS,  DLiPC,  DLiPE,   DMPA,   DMPC,   DMPE,   DMPG, 
    DMPS,   DNPA,   DNPC,   DNPE,   DNPG,   DNPS,   DOPA, 
    DOPC,   DOPE,   DOPG,  DOPP1,  DOPP2,  DOPP3,   DOPS, 
    DPPA,   DPPC,   DPPE,   DPPG,   DPPS,   DSPA,   DSPC, 
    DSPE,   DSPG,   DSPS,   DTPA,   DUPC,   DXPC,   DXPE, 
    DYPA,   DYPG,   DYPS,    ERG,   LLPA,   LLPC,   LLPE, 
    LLPS,   LPPC,    LSM,    NSM,    OSM,  PDOPC,  PDOPE, 
    PLPA,   PLPC,   PLPE,   PLPG,   PLPS,   POPA,   POPC, 
    POPE,   POPG,  POPP1,  POPP2,  POPP3,   POPS,    PSM, 
    SAPA,   SAPC,   SAPE,   SAPG,   SAPS,   SDPA,   SDPC, 
    SDPE,   SDPG,   SDPS,   SITO,   SLPA,   SLPC,   SLPE, 
    SLPG,   SLPS,   SOPA,   SOPC,   SOPE,   SOPG,   SOPS, 
    SSM,   STIG,   TIPA,   TSPC
    ---------------------------------------------------------------------------


This shows that there are two PDB collections in the built-in repository: one for water and ions and another for lipids.  The water/ion collection was created manually; these are very simple PDB files.  Any of those lipid residues can be referred to in a ``make_membrane_system`` task (See :ref:`subs_runtasks_make_membrane_system` and :ref:`example mper-tm symmetric bilayer` and :ref:`example mper-tm viral bilayer`).

The lipid collection was created using ``make-pdb-collection`` in the following way:  

.. code-block:: bash

    $ pestifer make-pdb-collection --streamID lipid
    $ pestifer make-pdb-collection --streamID lipid --substreamID cholesterol
    $ pestifer make-pdb-collection --streamID lipid --substreamID cholesterol --resname CHM1 --take-ic-from CHL1
    $ pestifer make-pdb-collection --streamID lipid --substreamID sphingo
    $ pestifer make-pdb-collection --streamID lipid --substreamID miscellaneous
    $ pestifer make-pdb-collection --streamID lipid --substreamID detergent --residueID C6DHPC
    $ tar zcf lipid.tgz lipid

The tarball ``lipid.tgz`` is the compressed PDB collection that pestifer uses, and it is contained in the ``resources`` data directory of the project.  The residue CHM1 does not have valid internal coordinates (ICs) because it is just a truncated version of the cholesterol residue CHL1, so we use the ``--take-ic-from`` option to copy the ICs from CHL1 to CHM1.

(As instructed in the current CHARMM force field, we use "model 1" for cholesterol.)

Contents of one lipid RESI entry in a PDB collection
++++++++++++++++++++++++++++++++++++++++++++++++++++

Each RESI in a PDB collection is represented by a subdirectory named after the RESI, and that subdirectory contains a set of PDB files that represent different configurations of the molecule.  For example, the ``DOPC`` RESI in the lipid collection has the following contents:

.. code-block:: text

    DOPC
    ├── DOPC-00.pdb
    ├── DOPC-01.pdb
    ├── DOPC-02.pdb
    ├── DOPC-03.pdb
    ├── DOPC-04.pdb
    ├── DOPC-05.pdb
    ├── DOPC-06.pdb
    ├── DOPC-07.pdb
    ├── DOPC-08.pdb
    ├── DOPC-09.pdb
    ├── DOPC-init.pdb
    ├── DOPC-init.psf
    ├── info.yaml
    └── init.tcl


The pdb files ``DOPC-00.pdb`` through ``DOPC-09.pdb`` are the 10 different configurations of the DOPC molecule.  The ``DOPC-init.pdb`` and ``DOPC-init.psf`` files are the initial coordinates and topology of the molecule, and the ``init.tcl`` file is a psfgen script used to generate those two files.  The ``info.yaml`` file contains metadata about the RESI, such as its long name and measurements of its dimensions that ``packmol`` needs:

.. code-block:: yaml

    charge: 0.0
    conformers:
    - head-tail-length: 27.324
      max-internal-length: 31.501
      pdb: POPC-00.pdb
    - head-tail-length: 28.827
      max-internal-length: 32.059
      pdb: POPC-01.pdb
    - head-tail-length: 28.67
      max-internal-length: 31.82
      pdb: POPC-02.pdb
    - head-tail-length: 28.222
      max-internal-length: 31.135
      pdb: POPC-03.pdb
    - head-tail-length: 27.051
      max-internal-length: 31.377
      pdb: POPC-04.pdb
    - head-tail-length: 26.786
      max-internal-length: 31.216
      pdb: POPC-05.pdb
    - head-tail-length: 27.72
      max-internal-length: 31.825
      pdb: POPC-06.pdb
    - head-tail-length: 27.918
      max-internal-length: 31.337
      pdb: POPC-07.pdb
    - head-tail-length: 27.752
      max-internal-length: 30.961
      pdb: POPC-08.pdb
    - head-tail-length: 27.942
      max-internal-length: 31.738
      pdb: POPC-09.pdb
    defined-in: top_all36_lipid.rtf
    parameters:
    - par_all36m_prot.prm
    - par_all36_na.prm
    - par_all36_cgenff.prm
    - toppar_all36_carb_glycopeptide.str
    - par_all36_carb.prm
    - toppar_water_ions.str
    - toppar_all36_prot_modify_res.str
    - par_all36_lipid.prm
    reference-atoms:
    heads:
    - name: N
      serial: 1
    tails:
    - name: C218
      serial: 88
    - name: C316
      serial: 131
    synonym: 3-palmitoyl-2-oleoyl-D-glycero-1-Phosphatidylcholine

Building your own PDB collections
+++++++++++++++++++++++++++++++++

Suppose you want to use lipid residues defined in the CHARMMFF stream file ``toppar_all36_lipid_yeast.str``; that is, you want PDBs for all the RESI's in the ``yeast`` substream. These are currently not part of the default PDB collection that comes with pestifer.  Consider the following commands:

.. code-block:: bash

    $ mkdir ~/my_pestifer_project
    $ cd ~/my_pestifer_project
    $ pestifer make-pdb-collection --streamID lipid --substreamID yeast --output-dir lipid-yeast

This will generated a directory ``~/my_pestifer_project/lipid-yeast/`` that contains the new PDB collection.  Each RESI subdirectory will contain 10 PDB files, each of which represents a different configuration of the molecule, along with an ``info.yaml`` file that contains important metadata about the RESI:  

.. code-block:: text

    lipid-yeast/
    ├── DYPC
    ├── DYPE
    ├── PYPE
    ├── YOPA
    ├── YOPC
    ├── YOPE
    └── YOPS

Each of these subdirectories contains the PDB files and metadata for that RESI.  For example, the ``DYPC`` subdirectory contains:

.. code-block:: text

    DYPC/
    ├── DYPC-00.pdb
    ├── DYPC-01.pdb
    ├── DYPC-02.pdb
    ├── DYPC-03.pdb
    ├── DYPC-04.pdb
    ├── DYPC-05.pdb
    ├── DYPC-06.pdb
    ├── DYPC-07.pdb
    ├── DYPC-08.pdb
    ├── DYPC-09.pdb
    ├── DYPC-init.pdb
    ├── DYPC-init.psf
    ├── info.yaml
    └── init.tcl

Suppose you want to use the PDB collection you just created in a ``make_membrane_system`` task.  You would need include the path in the ``pdbcollections`` list under the toplevel ``charmmff`` section:

.. code-block:: yaml

    charmmff:
      pdbcollections:
        - ~/my_pestifer_project/lipid-yeast
        