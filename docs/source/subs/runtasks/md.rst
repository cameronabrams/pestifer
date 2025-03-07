.. _subs_runtasks_md:

md 
--

An ``md`` task performs a NAMD simulation on the current PSF/PDB files; for example, those produced by the immediately previous ``psfgen`` task.

In most cases, it is only necessary to specify which ensemble (minimize, NVT, or NPT) you want to run in, and how many MD steps:

.. code-block:: yaml
   
   md:
     ensemble: NVT
     nsteps: 1000

The above ``md`` directive will run an NVT MD simulation at the default temperature of 300 K for 1000 time-steps.

Default NAMD config file values can be overwritten in the top-level :ref:`config_ref namd` directive of the pestifer config file.

A typical sequence of tasks for a new build might be a ``psfgen`` task, followed by an ``md`` task to perform a vacuum energy minimization, followed then by a short ``md`` task to thermalize the protein:

.. code-block:: yaml

   tasks:
     - psfgen:
         source:
          id: 6pti
     - md:
         ensemble: minimize
     - md:
         ensemble: NVT

If you are using GPU-resident namd3, I recommend setting the ``cpu-override`` flag to ``True`` in any ``md`` task run in vacuum.