.. _subs_runtasks_md:

md 
--

An ``md`` task performs a NAMD simulation on the current PSF/PDB files; for example, those produced by the immediately previous ``psfgen`` task.

In most cases, it is only necessary to specify which ensemble (minimize, NVT, NPT, or NPAT) you want to run in, and how many MD steps:

.. code-block:: yaml
   
   tasks:
     - ... (prior tasks here)
     - md:
         ensemble: NVT
         nsteps: 1000
     - ... (subsequent tasks here)

The above ``md`` directive will run an NVT MD simulation at the default temperature of 300 K for 1000 time-steps.  By default, any ensemble in which temperature is controlled will use a Langevin thermostat with a damping coefficient of 5/ps.  The default NPT and NPAT ensembles use a Langevin piston barostat with a target pressure of 1 atm and a piston period of 100 fs.  NPAT ensembles further restrict the ratio of the x and y dimensions of the box to maintain a constant ratio, which is useful for membrane simulations.

Default NAMD config file values can be overwritten in the top-level :ref:`config_ref namd` directive of the pestifer config file. The thermostat parameters can be overridden by specifying :ref:`config_ref namd thermostat` parameters, and the barostat parameters can be overridden by specifying :ref:`config_ref namd barostat` parameters.  The NPAT parameters can be overridden by specifying :ref:`config_ref namd membrane` parameters.

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