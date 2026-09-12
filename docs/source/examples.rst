.. _examples:

Examples
========

.. toctree::
   :hidden:
   :maxdepth: 1

   examples/01/bpti1
   examples/02/bpti2
   examples/03/bpti3
   examples/04/bpti4
   examples/05/hiv-protease
   examples/06/green-mamba-toxin
   examples/13/insulin-hexamer
   examples/19/sperm-whale-myoglobin
   examples/20/methylmalonyl-coa-mutase
   examples/21/groel-groes-adp
   examples/22/ferredoxin-fad
   examples/28/phosphoubiquitin
   examples/29/myristoylated-matrix
   examples/07/hiv-sosip-env-ectodomain1
   examples/08/hiv-sosip-env-ectodomain2
   examples/09/hiv-ad8-env-ectodomain
   examples/10/hiv-ae2-env-ectodomain
   examples/11/hiv-sosip-env-ectodomain3
   examples/12/hiv-sosip-env-ectodomain4
   examples/14/insulin-receptor-ectodomain
   examples/15/sars-cov2-S-BA2
   examples/16/hiv-mpertm3-membrane1
   examples/17/hiv-mpertm3-membrane2
   examples/23/subtilisin-dmso
   examples/24/subtilisin-acetone
   examples/26/subtilisin-acetonitrile
   examples/18/ecoli-polymerase
   examples/25/hiv-env-cd4-17b-liganded
   examples/27/ubiquitin-gfp-fusion

Fun with BPTI
-------------

Four builds from one source structure (`PDB 6pti <https://www.rcsb.org/structure/6PTI>`_,
bovine pancreatic trypsin inhibitor), each adding one more modification to a standard
fetch-psfgen-solvate-equilibrate workflow.  Start here.

- :doc:`Example 1 <examples/01/bpti1>` — baseline solvated build; the simplest complete pestifer workflow
- :doc:`Example 2 <examples/02/bpti2>` — heteroatom exclusion (phosphate ion); salty solvent; retaining crystal waters; ``validate`` task with custom parameters
- :doc:`Example 3 <examples/03/bpti3>` — point mutations (two shortcode formats) and disulfide bond deletion
- :doc:`Example 4 <examples/04/bpti4>` — introducing a new disulfide bond via mutations

Cofactors, metals and protonation
---------------------------------

Systems whose chemistry lives outside the standard amino-acid topologies: bound ligands,
heme and FAD cofactors, metal ions, oligomeric assemblies, post-translational modifications
carried through from the deposit, and residue protonation states assigned from the structure
rather than assumed.

- :doc:`Example 5 <examples/05/hiv-protease>` — HIV-1 protease dimer (1f7a); including small-molecule acetate ligands; reversing engineered mutations
- :doc:`Example 6 <examples/06/green-mamba-toxin>` — fasciculin 1 from green mamba snake venom (1fas); automatic ionization state assignment via ``pdb2pqr``
- :doc:`Example 13 <examples/13/insulin-hexamer>` — hexameric insulin (2ins); homo-oligomeric assembly from the asymmetric unit
- :doc:`Example 19 <examples/19/sperm-whale-myoglobin>` — sperm whale myoglobin (1mob); heme cofactor handled via standard CHARMM parameters
- :doc:`Example 20 <examples/20/methylmalonyl-coa-mutase>` — mitochondrial methylmalonyl-CoA mutase; AlphaFold model (UniProt P22033) as the input source
- :doc:`Example 21 <examples/21/groel-groes-adp>` — asymmetric GroEL/GroES chaperonin complex (1aon); 21-chain assembly with ADP ligands
- :doc:`Example 22 <examples/22/ferredoxin-fad>` — ferredoxin-NADP(H) reductase (2bgj); FAD cofactor; chain exclusion
- :doc:`Example 28 <examples/28/phosphoubiquitin>` — Ser65-phosphorylated ubiquitin (1ubq); installing a post-translational modification that is *not* in the input, by mutating the serine to ``SEP``
- :doc:`Example 29 <examples/29/myristoylated-matrix>` — myristoylated HIV-1 matrix protein (1uph); a lipid deposited as a separate ligand is fused into the residue CHARMM defines it as part of, keeping its NMR conformation

Model building: loops, glycans and cleavage
-------------------------------------------

Structures that arrive incomplete.  These examples build in unresolved loops, graft N-glycans,
substitute stub sequences across gaps too large to model, reverse engineered mutations, and
cleave a chain -- the work of turning a deposited structure into a simulatable molecule.

- :doc:`Example 7 <examples/07/hiv-sosip-env-ectodomain1>` — BG505 SOSIP 4zmj; glycan chainID reassignment; missing loop modeling; reversing SOSIP mutations
- :doc:`Example 8 <examples/08/hiv-sosip-env-ectodomain2>` — BG505 4tvp with liganded Fabs removed; multi-chain exclusion of Fab chains, waters, and ions
- :doc:`Example 9 <examples/09/hiv-ad8-env-ectodomain>` — 8fad; including a user-parameterized small molecule (CGenFF drug fragment)
- :doc:`Example 10 <examples/10/hiv-ae2-env-ectodomain>` — 8fae; same CGenFF-parameterized drug molecule, different structure
- :doc:`Example 11 <examples/11/hiv-sosip-env-ectodomain3>` — 7txd; loop substitutions; zero-occupancy C-terminus handling; Fab and sCD4 exclusion
- :doc:`Example 12 <examples/12/hiv-sosip-env-ectodomain4>` — 5vn3; Gly\ :sub:`3` stub substitutions for missing V1/V2 loops; sCD4 and Fab exclusion
- :doc:`Example 14 <examples/14/insulin-receptor-ectodomain>` — insulin receptor ectodomain (4zxb); large multi-chain complex with bound Fabs removed
- :doc:`Example 15 <examples/15/sars-cov2-S-BA2>` — BA.2 spike 7xix; N-glycan grafting from three donor PDB structures; loop modeling; furin cleavage

Membrane systems
----------------

Embedding a protein in a lipid bilayer, from a single-component patch to a multi-lipid
viral-mimetic membrane with an asymmetric leaflet composition.

- :doc:`Example 16 <examples/16/hiv-mpertm3-membrane1>` — gp41 MPER-TM trimer embedded in a pure DMPC bilayer
- :doc:`Example 17 <examples/17/hiv-mpertm3-membrane2>` — same trimer embedded in a multi-lipid viral-mimetic bilayer

Non-aqueous solvents
--------------------

The same enzyme (subtilisin Carlsberg) in three organic solvents, covering both a shipped
pre-equilibrated box and on-demand generation of one that is not shipped.

- :doc:`Example 23 <examples/23/subtilisin-dmso>` — subtilisin Carlsberg (1scd) in DMSO; non-aqueous solvation using a shipped pre-equilibrated solvent box
- :doc:`Example 24 <examples/24/subtilisin-acetone>` — subtilisin Carlsberg (1scd) in acetone; on-demand generation of a non-shipped CGenFF solvent box
- :doc:`Example 26 <examples/26/subtilisin-acetonitrile>` — subtilisin Carlsberg (1scd) in acetonitrile; another on-demand CGenFF solvent box

Composite and specialized workflows
-----------------------------------

Builds that compose several capabilities into one pipeline, including multi-script builds
in which helper configurations prepare inputs the main script consumes.

- :doc:`Example 18 <examples/18/ecoli-polymerase>` — *E. coli* replicative DNA polymerase (5fkw); DNA-protein complex with multiple chains and missing loops
- :doc:`Example 25 <examples/25/hiv-env-cd4-17b-liganded>` — sCD4/17b-liganded HIV-1 Env trimer (5vn3); completing partially-resolved 17b Fabs from a second structure (1gc1) with ``align`` + ``transfer_coords``, preserving the interface, then ``merge``
- :doc:`Example 27 <examples/27/ubiquitin-gfp-fusion>` — a ubiquitin-GFP fusion construct; the ``Cfusions`` modification appends GFP (1ema) to the C-terminus of ubiquitin (1ubq), auto-fetching and orienting the donor, building the GFP chromophore ``CRO`` in-chain, and aliasing selenomethionine
