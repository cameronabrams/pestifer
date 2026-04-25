.. _example ubiquitin:

Example 23: Human ubiquitin — alignment to a complex binding pose
-----------------------------------------------------------------

Human ubiquitin (`PDB ID 1ubq <https://www.rcsb.org/structure/1ubq>`_) is a
76-residue protein involved in targeting misfolded or damaged proteins for
proteasomal degradation.  It is one of the most conserved proteins across
eukaryotes and has been extensively studied by both NMR and X-ray
crystallography.

This example demonstrates the :ref:`align <subs_buildtasks_manipulate_mods_align>`
coordinate modification.  After the protein is built with ``psfgen`` and
subjected to a brief vacuum energy minimization, a ``manipulate`` task
rigidly superimposes the pipeline system onto a **different** crystal
structure — `PDB ID 1yd8 <https://www.rcsb.org/structure/1yd8>`_, a complex
of ubiquitin (chain A) bound to the GGA3 GAT domain (chain B).

The GGA proteins are monomeric clathrin adaptors that recognize
ubiquitinated cargo at the trans-Golgi network and direct it toward
late endosomes.  The GGA3 GAT domain contacts ubiquitin through a conserved
hydrophobic patch centred on Ile44, reorienting ubiquitin relative to its
free-solution crystal form.  Aligning the built ubiquitin system to chain A
of 1yd8 therefore places the protein in the orientation it adopts when
engaged with GGA3, which is the desired starting point if the goal is
to simulate ubiquitin in that recognition pose.

.. literalinclude:: ../../../../pestifer/resources/examples/ex23/inputs/ubiquitin.yaml
    :language: yaml

.. task-table:: ../../../../pestifer/resources/examples/ex23/inputs/ubiquitin.yaml


This build can be performed in a clean directory using:

.. code-block:: bash

   $ pestifer build-example 23

The key step is the ``manipulate`` task immediately after vacuum minimization:

.. code-block:: yaml

   - manipulate:
       mods:
         align:
           - ref_pdb: 1yd8.pdb
             mobile_sel: backbone
             ref_sel: "chain A and backbone"

``ref_pdb: 1yd8.pdb`` points to the 1yd8 crystal structure fetched earlier in
the build.  ``mobile_sel: backbone`` selects the backbone heavy atoms of the
built (and minimized) ubiquitin as the mobile set.  ``ref_sel: "chain A and
backbone"`` selects the corresponding atoms from the ubiquitin chain in 1yd8.
Both selections span the same 76 residues, so the atom-count congruency check
passes automatically.

After the fit matrix is computed, the **entire** pipeline system is rigidly
moved into the 1yd8 ubiquitin frame.  Solvation then proceeds in that frame,
ensuring the periodic box is oriented consistently with the GGA3-complex
crystal structure.

This pattern — build from one source, align to a different experimental
structure — is the intended use case for the ``align`` directive.  The
reference PDB can be any structure in a known or desired orientation: a
membrane-oriented model from the OPM database, a companion structure from a
different crystal form, or, as here, the binding pose from a protein–protein
complex.

.. raw:: html

    <div class="autogen-footer">
        <p>Example author: Cameron F. Abrams &nbsp;&nbsp;&nbsp; Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
    </div>
