.. _example myristoylated-matrix:

Myristoylated HIV-1 matrix protein
----------------------------------

`PDB ID 1uph <https://www.rcsb.org/structure/1UPH>`_ is an NMR structure of the HIV-1 matrix protein (MA, the N-terminal domain of Gag) carrying the myristoyl group that is covalently attached to its N-terminal glycine.  Myristoylation is what targets Gag to the plasma membrane; in this structure the fatty acyl chain is tucked against the protein rather than extended.  This example builds that modification *from the deposit* -- nothing is mutated -- and keeps the conformation the structure gives it.

The :ref:`phosphoubiquitin example <example phosphoubiquitin>` shows how to install a modification the input does not have.  This one is the other half: a modification the input already has, but written in a form CHARMM does not use.

**The deposit and CHARMM disagree about what a residue is.**  The PDB entry writes the myristate as a separate ligand, ``MYR A:1``, joined to Gly2 by a ``LINK`` record.  CHARMM defines N-myristoyl-glycine as a single residue, ``GLYM``, in ``toppar_all36_lipid_prot.str``.  pestifer reconciles the two when it reads the structure: a linked ``MYR`` is fused into its glycine as ``GLYM``, carrying the deposited heavy-atom coordinates and dropping the ligand's hydrogens for psfgen to rebuild.  The build log says so:

.. code-block:: text

    fused MYR A:1 into GLY A:2 as GLYM (15 heavy atoms carried, 27 hydrogens dropped for psfgen to rebuild)

Because the coordinates are carried rather than generated, the chain keeps its deposited conformation.  Measured on this build, model 1: all nine distal carbons (C6-C14) are packed against the protein, the chain stays within 0.45 Å mean (0.83 Å max) of the deposit after the first minimize, and the new amide bond is 1.35 Å against CHARMM's 1.345.  The fusion table also covers ``MYR`` on a lysine side chain (``LYSM``), though no example builds one.

Three things about this build are worth reading before adapting it.

**Pick one NMR model.**  ``source: model: 1`` selects the first conformer.  Without it, every model's coordinates are read, and the myristate arrives as twenty superimposed copies -- the log's "heavy atoms carried" count is 300 rather than 15.

**One parameter lives in an unrelated file.**  ``GLYM``'s topology is loaded automatically, but the angle, bond and dihedral terms at its acyl carbonyl (``C-CTL2``, ``C-CTL2-CTL2`` and four more) are defined only in ``toppar_all36_lipid_sphingo.str``, which nothing else about this system would load.  The config lists it:

.. code-block:: yaml

    charmmff:
      standard:
        str:
          - toppar_all36_lipid_sphingo.str

Entries under ``charmmff.standard.str`` are *added* to the default set, not substituted for it.  Leave the line out and the build stops before NAMD starts, naming the residue and every unresolved term:

.. code-block:: text

    ... needs force-field terms that none of this run's 12 parameter file(s) define:
      Unresolved bond term(s) [1]:
        'C-CTL2'  (e.g. residue GLYM A2)
      Unresolved angle term(s) [3]:
        'C-CTL2-CTL2'  (e.g. residue GLYM A2)
        ...

``LYSM`` and ``CYSL`` need the same stream, for the same six terms; the thioester ``CYSP`` and the prenylated ``CYSF`` and ``CYSG`` need nothing extra.  See :ref:`Post-translational modifications <post_translational_modifications>`.

**The chain starts without an N-terminal patch.**  ``GLYM`` acylates its own backbone nitrogen, so the usual ``NTER`` -- which turns that nitrogen into NH\ :sub:`3`\ :sup:`+` -- must not be applied.  pestifer writes ``first none`` for a segment that begins with ``GLYM``.  It has to: CHARMM's stream deliberately sets no terminal default, so the residue otherwise inherits whatever default the previously read topology file left in force, and with ``NTER`` in force psfgen builds an NH\ :sub:`3` nitrogen still bonded to the myristoyl carbonyl without complaint.

The ``validate`` task checks the outcome rather than assuming it: exactly one ``GLYM``, no free ``MYR`` left over, and 49 atoms in the fused residue (7 from glycine, 42 from the myristoyl chain).

.. literalinclude:: ../../../../pestifer/resources/examples/29/inputs/myristoylated-matrix.yaml
    :language: yaml

.. task-table:: ../../../../pestifer/resources/examples/29/inputs/myristoylated-matrix.yaml

.. raw:: html

    <div class="autogen-footer">
        <p>Example author: Cameron F. Abrams &nbsp;&nbsp;&nbsp; Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
    </div>
