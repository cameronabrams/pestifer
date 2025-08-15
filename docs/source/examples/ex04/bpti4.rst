.. _example bpti4:

Example 4: BPTI with a Mutated-in Disulfide Bond
------------------------------------------------

Using the ``mods`` subdirective, one can introduce new disulfides into an existing structure.  This example introduces a disulfide linking residues 11 and 34.

.. literalinclude:: ../../../../pestifer/resources/examples/ex04/inputs/bpti4.yaml
    :language: yaml

Note that this required first mutating the residues at positions 11 and 34 to cysteines, and *then* introducing the disulfide mod.

.. raw:: html

        <div class="autogen-footer">
            <p>Example author: Cameron F. Abrams&nbsp;&nbsp;&nbsp;Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
        </div>