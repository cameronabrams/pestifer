.. _subs_runtasks_continuation:

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
