.. _example ubiquitin:

Example 23: Human ubiquitin — crystal-frame alignment
------------------------------------------------------

Human ubiquitin (`PDB ID 1ubq <https://www.rcsb.org/structure/1ubq>`_) is a
76-residue protein involved in targeting misfolded or damaged proteins for
proteasomal degradation.  It is one of the most conserved proteins across
eukaryotes and has been extensively studied by both NMR and X-ray
crystallography.

This example demonstrates the :ref:`align <subs_buildtasks_manipulate_mods_align>`
coordinate modification.  After the protein is built with ``psfgen`` and
subjected to a brief vacuum energy minimization, a ``manipulate`` task
rigidly superimposes the pipeline system back onto the crystal structure
using a backbone least-RMSD fit, restoring the experimental coordinate frame
before solvation.

In practice, the reference PDB can be any structure in a known or desired
orientation — a previously equilibrated snapshot, a membrane-oriented model
from the OPM database, or a companion experimental structure from a different
crystal form.  Here, the reference is the crystal structure itself
(``1ubq.pdb``), which is written to the working directory by the ``fetch``
task and remains available throughout the build.

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
           - ref_pdb: 1ubq.pdb
             mobile_sel: backbone
             ref_sel: backbone

``ref_pdb: 1ubq.pdb`` points to the crystal-structure PDB file that was
fetched from the RCSB at the start of the build.  The ``backbone`` VMD
atom-selection string restricts the least-RMSD fit to the protein backbone
heavy atoms (N, Cα, C, O), excluding solvent.  Both the mobile and reference
selections must cover the same number of atoms; ``pestifer`` performs a
congruency check before invoking VMD's ``measure fit`` and will exit with an
error if they differ.

After the fit matrix is computed, the **entire** pipeline system (protein,
crystal waters, and any other molecules present) is rigidly moved into the
reference frame.  Solvation then proceeds in that frame, ensuring the final
periodic box is oriented consistently with the crystal structure.

.. raw:: html

    <div class="autogen-footer">
        <p>Example author: Cameron F. Abrams &nbsp;&nbsp;&nbsp; Contact: <a href="mailto:cfa22@drexel.edu">cfa22@drexel.edu</a></p>
    </div>
