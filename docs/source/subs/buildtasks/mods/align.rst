.. _subs_buildtasks_manipulate_mods_align:

align
-----

An ``align`` directive rigidly aligns the pipeline system (or a selection of it) to
a reference coordinate file by computing the least-RMSD rotation and translation
between two atom sets, then applying that transformation.  Internally it uses
VMD's ``measure fit`` command.

Parameters
++++++++++

``ref_pdb`` *(required)*
   Path to the reference PDB file.  The file must be present in the working
   directory when the task runs.

``ref_psf`` *(optional)*
   Path to a reference PSF file.  When provided, both the PSF and PDB are
   loaded together in VMD so that the full topology is available for selection
   purposes.  When omitted the reference is loaded as a plain PDB, which is
   sufficient for coordinate-based fitting.

``mobile_sel`` *(optional, default:* ``"all"`` *)*
   A VMD atomselect string identifying which atoms of the **pipeline system**
   are used to compute the fit.  The atom count and order must match those
   selected by ``ref_sel`` in the reference.

``ref_sel`` *(optional, default: inherits* ``mobile_sel`` *)*
   A VMD atomselect string identifying which atoms of the **reference system**
   are used to compute the fit.  When omitted, ``mobile_sel`` is reused,
   which is correct whenever both systems share the same topology.

``apply_to`` *(optional, default:* ``"all"`` *)*
   A VMD atomselect string identifying which atoms of the pipeline system are
   actually **moved** by the fitted transformation.  This allows you to compute
   the fit on a subset (e.g. backbone only) while translating and rotating
   the entire system.

Congruency check
++++++++++++++++

Before computing the fit, the generated VMD script verifies that ``mobile_sel``
and ``ref_sel`` select the same number of atoms.  If they differ, VMD prints a
descriptive error message and exits with a non-zero code, which pestifer
surfaces as a task failure.

.. note::

   ``measure fit`` requires both selections to contain the same number of atoms
   **in the same order**.  Atom-count equality is checked automatically; atom
   ordering is the user's responsibility.  When in doubt, use consistent
   selection strings and load the reference with its PSF.

Examples
++++++++

Align the entire pipeline system to a reference PDB (PDB-only reference, fit
and move all atoms):

.. code-block:: yaml

   tasks:
     - manipulate:
         mods:
           align:
             - ref_pdb: reference.pdb

Fit on backbone atoms only, move the whole system, using a reference that
has its own topology:

.. code-block:: yaml

   tasks:
     - manipulate:
         mods:
           align:
             - ref_pdb: reference.pdb
               ref_psf: reference.psf
               mobile_sel: "backbone"
               apply_to: "all"

Fit protein chain A to a reference chain using different selection strings
for mobile and reference:

.. code-block:: yaml

   tasks:
     - manipulate:
         mods:
           align:
             - ref_pdb: reference.pdb
               mobile_sel: "chain A and backbone"
               ref_sel: "chain B and backbone"
               apply_to: "all"

Workflow note
+++++++++++++

The ``align`` directive moves the **pipeline system** to match the reference.
If instead you need to align an incoming external system to the pipeline
system (for example, before a :ref:`merge <subs_buildtasks_merge>`), that
external system must first have its own PSF built in a separate pestifer run,
because ``merge`` requires a PSF for every input system.  The alignment is
a coordinate-only operation so the PSF from that separate run can be reused
directly; there is no need to rebuild it after alignment.
