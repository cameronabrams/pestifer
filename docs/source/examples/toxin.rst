.. _example toxin:

Example 6: Green Mamba Toxin at pH 7.0
--------------------------------------

You can use a ``pdb2pqr`` task to calculate and assign protonation states to ionizable residues. This example uses the green mamba toxin structure from `PDB ID 1fas <https://www.rcsb.org/structure/1fas>`_ at pH 7.0.

.. literalinclude:: ../../../pestifer/resources/examples/toxin.yaml
    :language: yaml

Note that the ``pdb2pqr`` task is invoked immediately after the ``psfgen`` task but *before* the ``solvate`` task. The ``pdb2pqr`` task is documented at :ref:`subs_runtasks_pdb2pqr`.
