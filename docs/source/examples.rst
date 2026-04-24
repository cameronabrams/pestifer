.. _examples:

Examples
========

.. toctree::
   :hidden:
   :maxdepth: 1

   examples/ex01/bpti1
   examples/ex02/bpti2
   examples/ex03/bpti3
   examples/ex04/bpti4
   examples/ex05/hiv-protease
   examples/ex06/green-mamba-toxin
   examples/ex07/hiv-sosip-env-ectodomain1
   examples/ex08/hiv-sosip-env-ectodomain2
   examples/ex09/hiv-ad8-env-ectodomain
   examples/ex10/hiv-ae2-env-ectodomain
   examples/ex11/hiv-sosip-env-ectodomain3
   examples/ex12/hiv-sosip-env-ectodomain4
   examples/ex13/insulin-hexamer
   examples/ex14/insulin-receptor-ectodomain
   examples/ex15/sars-cov2-S-BA2
   examples/ex16/hiv-mpertm3-membrane1
   examples/ex17/hiv-mpertm3-membrane2
   examples/ex18/ecoli-polymerase
   examples/ex19/sperm-whale-myoglobin
   examples/ex20/methylmalonyl-coa-mutase
   examples/ex21/groel-groes-adp
   examples/ex22/ferredoxin-fad
   examples/ex23/ubiquitin

BPTI Series (Examples 1–4)
---------------------------

All four examples build from the same source structure (`PDB 6pti
<https://www.rcsb.org/structure/6PTI>`_, bovine pancreatic trypsin inhibitor),
demonstrating progressively more complex modifications to a standard
fetch-psfgen-solvate-equilibrate workflow.

- :doc:`Example 1 <examples/ex01/bpti1>` — baseline solvated build; the simplest complete pestifer workflow
- :doc:`Example 2 <examples/ex02/bpti2>` — heteroatom exclusion (phosphate ion); salty solvent; retaining crystal waters; ``validate`` task with custom parameters
- :doc:`Example 3 <examples/ex03/bpti3>` — point mutations (two shortcode formats) and disulfide bond deletion
- :doc:`Example 4 <examples/ex04/bpti4>` — introducing a new disulfide bond via mutations

Small Proteins (Examples 5–6)
------------------------------

Straightforward builds of small, single-chain or homodimeric proteins,
illustrating ligand inclusion and the reversal of engineered mutations.

- :doc:`Example 5 <examples/ex05/hiv-protease>` — HIV-1 protease dimer (1f7a); including small-molecule acetate ligands; reversing engineered mutations
- :doc:`Example 6 <examples/ex06/green-mamba-toxin>` — fasciculin 1 from green mamba snake venom (1fas); automatic ionization state assignment via ``pdb2pqr``

HIV-1 Env Ectodomains (Examples 7–12)
--------------------------------------

A family of HIV-1 Env SOSIP trimer builds drawn from several crystal
and cryo-EM structures, collectively demonstrating glycan chain management, Fab and
ligand exclusion, missing-loop modeling, zero-occupancy terminus handling,
user-parameterized small molecules, and stub substitutions for absent loops.

- :doc:`Example 7 <examples/ex07/hiv-sosip-env-ectodomain1>` — BG505 SOSIP 4zmj; glycan chainID reassignment; missing loop modeling; reversing SOSIP mutations
- :doc:`Example 8 <examples/ex08/hiv-sosip-env-ectodomain2>` — BG505 4tvp with liganded Fabs removed; multi-chain exclusion of Fab chains, waters, and ions
- :doc:`Example 9 <examples/ex09/hiv-ad8-env-ectodomain>` — 8fad; including a user-parameterized small molecule (CGenFF drug fragment)
- :doc:`Example 10 <examples/ex10/hiv-ae2-env-ectodomain>` — 8fae; same CGenFF-parameterized drug molecule, different structure
- :doc:`Example 11 <examples/ex11/hiv-sosip-env-ectodomain3>` — 7txd; loop substitutions; zero-occupancy C-terminus handling; Fab and sCD4 exclusion
- :doc:`Example 12 <examples/ex12/hiv-sosip-env-ectodomain4>` — 5vn3; Gly\ :sub:`3` stub substitutions for missing V1/V2 loops; sCD4 and Fab exclusion

Insulin Systems (Examples 13–14)
---------------------------------

Builds involving insulin and its receptor, illustrating homo-oligomeric
assembly from an asymmetric unit and Fab exclusion from a large multi-chain
complex.

- :doc:`Example 13 <examples/ex13/insulin-hexamer>` — hexameric insulin (2ins); homo-oligomeric assembly from the asymmetric unit
- :doc:`Example 14 <examples/ex14/insulin-receptor-ectodomain>` — insulin receptor ectodomain (4zxb); large multi-chain complex with bound Fabs removed

SARS-CoV-2 Spike (Example 15)
------------------------------

A fully glycosylated BA.2 SARS-CoV-2 spike trimer build, demonstrating
glycan grafting from donor structures, missing-loop modeling, and furin
cleavage-site handling.

- :doc:`Example 15 <examples/ex15/sars-cov2-S-BA2>` — BA.2 spike 7xix; N-glycan grafting from three donor PDB structures; loop modeling; furin cleavage

Membrane-Embedded Systems (Examples 16–17)
-------------------------------------------

Both examples embed the same HIV-1 gp41 MPER-TM trimer (6e8w) in a lipid
bilayer, contrasting a simple single-lipid membrane with a multi-component
viral-mimetic bilayer.

- :doc:`Example 16 <examples/ex16/hiv-mpertm3-membrane1>` — gp41 MPER-TM trimer embedded in a pure DMPC bilayer
- :doc:`Example 17 <examples/ex17/hiv-mpertm3-membrane2>` — same trimer embedded in a multi-lipid viral-mimetic bilayer

Diverse Applications (Examples 18–22)
--------------------------------------

A collection of builds that highlight additional pestifer capabilities:
AlphaFold models as input, cofactor-containing enzymes, DNA-protein complexes,
and very large multi-chain assemblies.

- :doc:`Example 18 <examples/ex18/ecoli-polymerase>` — *E. coli* replicative DNA polymerase (5fkw); DNA-protein complex with multiple chains and missing loops
- :doc:`Example 19 <examples/ex19/sperm-whale-myoglobin>` — sperm whale myoglobin (1mob); heme cofactor handled via standard CHARMM parameters
- :doc:`Example 20 <examples/ex20/methylmalonyl-coa-mutase>` — mitochondrial methylmalonyl-CoA mutase; AlphaFold model (UniProt P22033) as the input source
- :doc:`Example 21 <examples/ex21/groel-groes-adp>` — asymmetric GroEL/GroES chaperonin complex (1aon); 21-chain assembly with ADP ligands
- :doc:`Example 22 <examples/ex22/ferredoxin-fad>` — ferredoxin-NADP(H) reductase (2bgj); FAD cofactor; chain exclusion

Coordinate Modifications (Example 23)
--------------------------------------

Builds that showcase the :ref:`manipulate <subs_buildtasks_manipulate>` task's
coordinate-only modification capabilities, which reposition atoms after the PSF
has been built without altering system topology.

- :doc:`Example 23 <examples/ex23/ubiquitin>` — human ubiquitin (1ubq); post-minimization backbone alignment to the crystal coordinate frame using the ``align`` mod
