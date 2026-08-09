.. _subs_buildtasks_continuation:

continuation
------------

A ``continuation`` task allows you begin the run from an already prepared system state, represented by PDB, PSF, and XSC files.

For example, to begin a run from a state represented by ``my_6pti.psf``, ``my_6pti.pdb``, and ``my_6pti.xsc``, your first
task would look like this:

.. code-block:: yaml
        
    tasks:
    - continuation:
        psf: my_6pti.psf
        pdb: my_6pti.pdb
        xsc: my_6pti.xsc

A ``psf`` and either a ``pdb`` or a ``coor`` are required; ``xsc`` (periodic box), ``coor`` (binary coordinates), and ``vel`` (binary velocities) are optional.  All files provided should correspond to the same state.

Importing a foreign PSF
=======================

The incoming state does not have to be one pestifer built.  A system prepared by CHARMM-GUI, by another tool, or by an older pestifer version enters the pipeline through ``continuation`` and flows on into MD, membrane embedding, packaging -- and a following :ref:`psfgen <subs_buildtasks_psfgen>` task can *edit* it.

Because a foreign PSF may have been built against a different CHARMM release than the one pestifer is using, ``continuation`` runs a **force-field consistency preflight** on it, controlled by ``verify_parameters`` (default ``true``):

.. code-block:: yaml

    tasks:
    - continuation:
        psf: from_charmm_gui.psf
        pdb: from_charmm_gui.pdb
        xsc: from_charmm_gui.xsc
        verify_parameters: true

The preflight parses the PSF once and does two things.  First, it **inventories** the structure: it classifies each segment's residues and reports a per-segid composition summary, warning about any residue name it cannot map to a known segment type.  Second, it **verifies force-field coverage**: every atom type, and every bond, angle, dihedral, and improper type-tuple in the PSF, must resolve against the build's CHARMM release.  If any does not, the build stops immediately with a message naming the residue and the offending term.

The point is to fail *early and legibly*.  Without the preflight, a missing parameter surfaces as a cryptic ``DIDN'T FIND vdW PARAMETER`` abort from NAMD, potentially hours into a build.

.. note::

   Pestifer's own internal continuations (those it generates between stages of a multi-pass build) set ``verify_parameters: false``, since the PSF was just built against the same release; this keeps the fast, zero-parse path.  Set it to ``false`` yourself only if you are confident the PSF matches your force field and want to skip the parse.

If the incoming PSF needs stream files (a glycan or lipid topology stream, say), pestifer resolves them from the ``REMARKS topology`` records the PSF carries, copying them into the run directory so a downstream ``psfgen`` can read them.

.. note::

   If any of the specified files reside outside the current working directory, Pestifer will
   automatically copy them into the CWD before registering them.  This ensures that subsequent
   tasks — which write their outputs alongside the current state files — always operate within
   a single directory, which is important for multi-pass workflows where the input state was
   produced in a different directory.
