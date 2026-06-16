.. _subs_buildtasks_mdtask:

md
--

An ``mdtask`` performs a NAMD simulation on the current state files; for example, those produced by the immediately previous ``psfgen`` task.

In most cases, it is only necessary to specify which ensemble (minimize, NVT, NPT, or NPgT) you want to run in, and how many MD steps:

.. code-block:: yaml
   
   tasks:
     - ... (prior tasks here)
     - md:
         ensemble: NVT
         nsteps: 1000
     - ... (subsequent tasks here)

The above ``md`` directive will run an NVT MD simulation at the default temperature of 300 K for 1000 time-steps.  By default, any ensemble in which temperature is controlled will use a Langevin thermostat with a damping coefficient of 5/ps.  The default NPT and NPgT ensembles use a Langevin piston barostat with a target pressure of 1 atm.  The **NPgT** ensemble further constrains the x and y box dimensions to a **constant ratio** (they scale together while the cross-sectional area is free to change, ``useConstantRatio``), which is useful for membrane simulations; NPT instead scales the cell isotropically.

.. note::

   ``NPgT`` is the preferred name for this constant-x:y-ratio ensemble.  Pestifer
   historically called it ``NPAT``, but that label has always set
   ``useFlexibleCell`` + ``useConstantRatio`` (constant *ratio*), **not** the
   ``useConstantArea`` that the "A = area" in the NPAT acronym implies.  ``NPAT`` is
   still accepted (and behaves identically to ``NPgT``) but now emits a deprecation
   warning: a future version will reserve ``NPAT`` for true constant-area runs.  The
   names ``NPgT`` (preferred), ``NPGT``, and ``npgt`` are all equivalent.

Default NAMD config file values can be overwritten in the top-level :ref:`config_ref namd` directive of the pestifer config file. The thermostat parameters can be overridden by specifying :ref:`config_ref namd thermostat` parameters, and the barostat parameters can be overridden by specifying :ref:`config_ref namd barostat` parameters.  The NPgT parameters can be overridden by specifying :ref:`config_ref namd membrane` parameters.

A typical sequence of tasks for a new build might be a ``fetch`` task, then a ``psfgen`` task, followed by an ``md`` task to perform a vacuum energy minimization, followed then by a short ``md`` task to thermalize the protein:

.. code-block:: yaml

   tasks:
     - fetch:
         sourceID: 6pti
     - psfgen:
     - md:
         ensemble: minimize
     - md:
         ensemble: NVT

If you are using GPU-resident namd3, I recommend setting the ``cpu-override`` flag to ``True`` in any ``md`` task run in vacuum.