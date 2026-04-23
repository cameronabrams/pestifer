.. _subs_buildtasks_manipulate_mods_transfer_coords:

transfer_coords
---------------

A ``transfer_coords`` directive copies atomic coordinates from a selected
subset of a **donor** PDB file directly onto the matching atoms of the
pipeline system, overwriting their positions without modifying any
topology.  This makes it suitable for non-rigid changes — re-modeled
loops, isolated energy minimization of a ligand or domain, manual
repositioning in an external editor — that cannot be represented as a
single rigid-body transformation.

An optional **pre-alignment** step can be requested.  When the donor was
prepared in a different coordinate frame (e.g. the loop was minimized as
a standalone fragment), supplying ``align_donor_sel`` and
``align_mobile_sel`` causes the entire donor to be rigidly fitted to the
pipeline system before the transfer, so that the transferred coordinates
land in the correct frame.

Parameters
++++++++++

``donor_pdb`` *(required)*
   Path to the donor PDB file.

``donor_psf`` *(optional)*
   Path to a donor PSF file.  When provided both files are loaded so that
   topology attributes (segment name, residue type, etc.) can be used in
   VMD atomselect strings.  When omitted the donor is loaded as a plain
   PDB.

``donor_sel`` *(required)*
   VMD atomselect string identifying which atoms to copy **from** the donor.

``mobile_sel`` *(required)*
   VMD atomselect string identifying which atoms of the pipeline system to
   **overwrite**.  Must select the same number of atoms as ``donor_sel``
   in the same order.

``align_donor_sel`` *(optional)*
   VMD atomselect string for donor atoms used in the pre-alignment fit.
   Must be supplied together with ``align_mobile_sel``; providing only one
   of the two is an error.

``align_mobile_sel`` *(optional)*
   VMD atomselect string for pipeline atoms used in the pre-alignment fit.
   Must be supplied together with ``align_donor_sel``.

.. note::

   Both the alignment pair and the transfer pair are checked for atom-count
   congruency before use.  If either check fails, VMD prints a descriptive
   error and exits with a non-zero code, which pestifer surfaces as a task
   failure.  Atom ordering within each congruent pair is the user's
   responsibility.

   Unlike :ref:`align <subs_buildtasks_manipulate_mods_align>`, which moves
   the **pipeline** system to match a reference, ``transfer_coords`` moves
   **donor** coordinates into the pipeline system.

Examples
++++++++

Copy coordinates of a re-modeled loop (same coordinate frame, no
pre-alignment needed):

.. code-block:: yaml

   tasks:
     - manipulate:
         mods:
           transfer_coords:
             - donor_pdb: remodeled_loop.pdb
               donor_sel: "chain A and resid 45 to 60"
               mobile_sel: "chain A and resid 45 to 60"

Inject coordinates from a donor that was prepared as a standalone
fragment (different coordinate frame — pre-alignment required):

.. code-block:: yaml

   tasks:
     - manipulate:
         mods:
           transfer_coords:
             - donor_pdb: minimized_loop.pdb
               donor_sel: "resid 45 to 60"
               mobile_sel: "chain A and resid 45 to 60"
               align_donor_sel: "backbone"
               align_mobile_sel: "chain A and backbone"

Combining with :ref:`align <subs_buildtasks_manipulate_mods_align>` in one
task (first align the whole system, then inject a refined loop):

.. code-block:: yaml

   tasks:
     - manipulate:
         mods:
           align:
             - ref_pdb: reference.pdb
               mobile_sel: "backbone"
           transfer_coords:
             - donor_pdb: refined_loop.pdb
               donor_sel: "resid 45 to 60"
               mobile_sel: "chain A and resid 45 to 60"
               align_donor_sel: "backbone"
               align_mobile_sel: "chain A and backbone"
