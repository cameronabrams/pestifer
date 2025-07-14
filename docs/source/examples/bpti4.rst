.. _example bpti4:

Example 4: BPTI with a Mutated-in Disulfide Bond
------------------------------------------------

Using the ``mods`` subdirective, one can introduce new disulfides into an existing structure.  This example introduces a disulfide linking residues 11 and 34.

.. literalinclude:: ../../../pestifer/resources/examples/bpti4/bpti4.yaml
    :language: yaml

Note that this required first mutating the residues at positions 11 and 34 to cysteines, and *then* introducing the disulfide mod.
