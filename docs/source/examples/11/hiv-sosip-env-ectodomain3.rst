.. _example env 7txd:

Example 11: Open, Symmetric D9/CD4-liganded HIV-1 SOSIP Env Ectodomain Trimer (ligands removed)
-----------------------------------------------------------------------------------------------

`PDB ID 7txd <https://www.rcsb.org/structure/7txd>`_ is a structure of the HIV-1 Env BG505 SOSIP-664 ectodomain trimer in an open conformation, with the two-domain soluble CD4 construct (sCD4) and D9 Fab fragment bound.  This structure is missing the so-called V1/2 loop between residues 100 and 162 on each gp120.  Here, we substitute that span of residues with a literal triglycine (Gly-Gly-Gly) sequence--this ensures psfgen doesn't just pop standard N- and C-terminal patches across that gap.  This is an alternative to what we did in :ref:`4zmj <example env 4zmj>`, where we opted to build in missing residues following sequence.   We also model in residues missing from gp41's and fold them into alpha-helices.  The ligands are removed from the structure, and the system is solvated in a water box.

Like the other SOSIP examples here, the configuration reverts the engineered stabilizing changes back to the wild-type sequence: the ``SOS`` disulfide (A501C and T605C, which staples gp120 to gp41) and the ``IP`` substitution (I559P).  This is not automatic -- it happens because the ``mutations`` block asks for it, and a build without that block is of the stabilized construct rather than of wild-type Env.  The names above are the conventional HXB2 ones; each configuration addresses the residues in its own structure's numbering, so the indices in the ``mutations`` block below will not always match them.

.. literalinclude:: ../../../../pestifer/resources/examples/11/inputs/hiv-sosip-env-ectodomain3.yaml
    :language: yaml

.. task-table:: ../../../../pestifer/resources/examples/11/inputs/hiv-sosip-env-ectodomain3.yaml


.. raw:: html

    <div class="autogen-footer">
        <p>Example author: Cameron F. Abrams &nbsp;&nbsp;&nbsp; Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
    </div>