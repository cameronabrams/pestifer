.. _example 4:

Example 4: BPTI with an Extra Disulfide Added
---------------------------------------------

Using the ``mods`` subdirective, one can introduce new disulfides into an existing structure.  This example introduces a disulfide linking residues 11 and 34.

.. literalinclude:: ../../../pestifer/resources/examples/04-bpti-mutant-extradisulf.yaml
    :language: yaml

Note that this required first mutating the residues at positions 11 and 34 to cysteines, and *then* introducing the disulfide mod.
