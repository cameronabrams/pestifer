.. _example bpti no md:

BPTI Built and Packaged Without an MD Step
-------------------------------------------

This example is here for its **shape**, not its science.  Every other bundled example runs at
least one ``md``, ``density_equilibrate`` or ``membrane_equilibrate`` task before it terminates.
This one runs none: it fetches a structure, builds the PSF, checks it, and packages the result.

That is a perfectly ordinary thing for a user to want -- a structure that needs no relaxation, or
a run whose only job is to produce a package -- and it takes a **different path** through
:ref:`terminate <subs_buildtasks_terminate>`.  The standard CHARMM parameter files are staged as
artifacts by MD tasks, so when a build runs none of them, ``terminate`` has to reach for the
standard set itself in order to write the consolidated ``*_minimal.prm``.  A real defect lived on
that branch undetected through v3.24.0 -- the packaged parameter file was merged in a *different
order* from the one a NAMD run would have used, so a duplicated force-field term could resolve
differently inside the tarball than it did during the build.

It went unnoticed because the example set is pestifer's only end-to-end coverage, and no example
had this shape.  That is the gap this example closes: not an untested *option*, but an untested
*task sequence*, which no amount of per-option coverage would have found.

.. literalinclude:: ../../../../pestifer/resources/examples/32/inputs/bpti-no-md.yaml
    :language: yaml

.. task-table:: ../../../../pestifer/resources/examples/32/inputs/bpti-no-md.yaml

The build takes about three seconds, which is the point: a shape worth covering should be cheap
enough that carrying it costs nothing.

Two options are exercised here that no other example turns on:

``psfgen.source.reserialize``
    Renumber atom serials from 1 when the PDB is written, rather than preserving the input's
    numbering.

``terminate.cleanup``
    Set to ``false``, so the intermediate files are left in place instead of being swept at the
    end.  Useful when you want to inspect what a build produced; the default ``true`` is what the
    other examples use.

What to check when it finishes
==============================

The package is the deliverable, and it should be self-contained -- every file named inside the
packaged NAMD config should be present in the tarball beside it:

.. code-block:: console

  $ tar tzf prod_6pti_nomd.tar.gz
  prod_6pti_nomd/my_6pti_nomd.psf
  prod_6pti_nomd/my_6pti_nomd.pdb
  prod_6pti_nomd/prod_6pti_nomd.namd
  prod_6pti_nomd/my_6pti_nomd_minimal.prm

and the consolidated parameter file should carry real parameters rather than an empty
``NONBONDED`` section -- this build produces 37 vdW entries, 87 bonds, 260 angles and 541
dihedrals.  An empty section is the symptom that the standard parameters were never picked up,
which is exactly what this path is at risk of.
