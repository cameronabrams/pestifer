.. _charmmff_customizations:

CHARMM force-field customizations
=================================

Pestifer bundles a CHARMM force field for each supported release under
``pestifer/resources/charmmff/<release>/`` (currently ``feb26``). The bundled
toppar tarball is **trimmed** to the files pestifer actually loads — the top-level
topology/parameter/stream files plus the ``stream/`` subtree — via
``scripts/trim_toppar.sh``; the upstream subdirectories pestifer never reads
(``metals/``, ``drude/``, ``non_charmm/``, implicit-solvent sets, etc.) are dropped.
Beyond the stock CHARMM toppar files, there are two kinds of local additions:

* ``<release>/patches/`` — **records of corrections** applied by hand to specific
  upstream toppar files. These live **per release**, since each patch targets a
  particular release's toppar file.
* ``custom/`` — **extra** topology / parameter / stream files that are not part of
  the base CHARMM release. These are release-independent, so they live in a **single
  shared** ``charmmff/custom/`` directory that is a sibling of the release directories
  (not duplicated under each one).

.. note::

   The ``.patch`` files under ``patches/`` are **provenance only** — nothing in
   pestifer applies them at build or package-regeneration time. The corresponding
   upstream ``.str`` file shipped in the release already contains the fix; the
   ``.patch`` is kept so the exact change (and the reason for it) is recoverable.

Corrections to upstream files (``patches/``)
--------------------------------------------

``toppar_all36_lipid_cardiolipin.patch``
   Corrects internal-coordinate (IC) table entries for cardiolipin acyl-chain
   hydrogens. For the second hydrogen of several methylene pairs (``H12A``/``H12B``
   and ``H13A``/``H13B`` on both the A- and C-chains), the released file repeats the
   *previous* IC line's leading atom instead of naming the partner hydrogen — e.g.
   ``IC CA13 CA11 *CA12 H12B`` where it should read ``IC H12A CA11 *CA12 H12B``. With
   the wrong leading atom, ``guesscoord`` builds the two hydrogens of the pair at the
   same position. The patch sets the correct leading atom so each hydrogen is placed
   distinctly.

``toppar_all36_prot_modify_res.patch``
   Corrects the coordinating nitrogen in the two histidine–zinc covalent-link patches.
   By name and intent, ``ZNHD`` bonds the zinc to a histidine's **delta-1 nitrogen**
   (``ND1``) and ``ZNHE`` to its **epsilon-2 nitrogen** (``NE2``). The released file
   has the two bonds swapped (``ZNHD → NE2``, ``ZNHE → ND1``); pestifer's copy restores
   them so each patch links the zinc through the nitrogen it is named for. These
   patches are used in the insulin-hexamer example build (Example 13), whose zinc ions
   are coordinated by histidines.

Added definitions (``custom/``)
-------------------------------

These files supply residues, patches, parameters, and ions not present in the base
CHARMM release. They live in the shared ``charmmff/custom/`` directory and are loaded
alongside the stock toppar for whichever release is active.

``pestifer.top``
   A small topology file of pestifer-specific patches:

   * a patch converting a standard up-puckered proline to a down-puckered proline
     (``PROD``);
   * several terminus-"undoing" patches (to reverse standard N-/C-terminal patching);
   * a deprotonated tyrosine.

``toppar_all36_moreions.str``
   Additional monatomic-ion definitions, taken from CGenFF.

``vcg-paramchem-ic.str``
   A CGenFF/ParamChem stream file for **Polysorbate-80** (``RESI VCG``), including the
   internal coordinates needed to build it.

``83G-cgenff.str``
   CGenFF/ParamChem parameters for ligand ``83G`` — the HIV-1 entry inhibitor
   **BMS-378806**, present in the HIV-1 Env ectodomain examples (see Examples 9 and 10).

``LF0-cgenff.str``
   CGenFF/ParamChem parameters for ligand ``LF0`` — an **HIV-1 integrase inhibitor**
   (see PDB entries ``4ID1``, ``4LH5``, ``4TSX``, ``8A1P``, ``8USY``, ``8V0Z``).

Contributing new custom definitions
-----------------------------------

The ``83G`` and ``LF0`` streams are examples of **built-in custom** residue
definitions. To add your own from any file containing a CHARMM ``RESI`` block (a
``.str``, ``.rtf``, or ``.top`` file), use ``modify-package`` — see
:ref:`subs modify-package` — which installs it into the shared ``charmmff/custom/``
directory and records the change in the modifications ledger.
