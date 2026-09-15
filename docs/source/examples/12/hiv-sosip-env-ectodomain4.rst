.. _example env 5vn3:

Open, Symmetric 17b/CD4-liganded HIV-1 B41 SOSIP Env Ectodomain Trimer
----------------------------------------------------------------------

`PDB ID: 5vn3 <https://www.rcsb.org/structure/5vn3>`_ represents one of the earliest structures of an open form of the HIV-1 Env ectodomain trimer.  Here we prepare a system with just the ectodomain, omitting the sCD4 and 17b chains.

Each gp120 is missing its V1/V2 loop, and rather than model that span in, we replace it with a literal triglycine (Gly-Gly-Gly) stub -- enough to stop psfgen capping the gap with standard N- and C-terminal patches, without inventing a loop conformation the structure does not support.  The :ref:`7txd example <example env 7txd>` does the same thing; the :ref:`4zmj example <example env 4zmj>` takes the other route and builds the missing residues in following sequence.

Like the other SOSIP examples here, the configuration reverts the engineered stabilizing changes back to the wild-type sequence: the ``SOS`` disulfide (A501C and T605C, which staples gp120 to gp41) and the ``IP`` substitution (I559P).  This is not automatic -- it happens because the ``mutations`` block asks for it, and a build without that block is of the stabilized construct rather than of wild-type Env.  The names above are the conventional HXB2 ones; each configuration addresses the residues in its own structure's numbering, so the indices in the ``mutations`` block below will not always match them.

.. literalinclude:: ../../../../pestifer/resources/examples/12/inputs/hiv-sosip-env-ectodomain4.yaml
    :language: yaml

.. task-table:: ../../../../pestifer/resources/examples/12/inputs/hiv-sosip-env-ectodomain4.yaml


.. raw:: html

    <div class="autogen-footer">
        <p>Example author: Cameron F. Abrams &nbsp;&nbsp;&nbsp; Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
    </div>