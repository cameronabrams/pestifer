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

**All three** of these files are required, and they should all correspond to the same state.

.. note::

   If any of the specified files reside outside the current working directory, Pestifer will
   automatically copy them into the CWD before registering them.  This ensures that subsequent
   tasks — which write their outputs alongside the current state files — always operate within
   a single directory, which is important for multi-pass workflows where the input state was
   produced in a different directory.
