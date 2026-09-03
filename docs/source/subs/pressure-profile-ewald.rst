.. _subs_pressure_profile_ewald:

pressure-profile-ewald
----------------------

The ``pressure-profile-ewald`` subcommand reconstructs a **complete** NAMD pressure profile — PME reciprocal-space term included — by replaying a saved trajectory twice and summing the two halves frame by frame.

Why two runs are needed
~~~~~~~~~~~~~~~~~~~~~~~

A run with ``pressureProfile on`` reports, in NAMD's own words, the "kinetic, bonded and nonbonded (but not reciprocal space) contributions".  The missing reciprocal term is not a small correction: on one replica of :ref:`example 16 <example mper-tm symmetric bilayer>` the real-space integral is :math:`-16.63` mN/m and the reciprocal part :math:`+26.50` mN/m, so including it *flips the sign* of the result.

The obvious fix — switching on ``pressureProfileEwald`` as well — does not work, because NAMD makes the two mutually exclusive.  From ``SimParameters.C`` (3.0.3, line 6699):

.. code-block:: c

   // if Ewald is on, only calculate Ewald
   if (pressureProfileEwaldOn)
     pressureProfileOn = 0;

So a run with ``pressureProfileEwald on`` reports the reciprocal contribution **alone**.  There is no single-run setting that yields the whole profile; it must be assembled from two.

Nor is it practical to compute the Ewald half during the production run.  NAMD registers ``ComputeEwald`` as an ordinary compute gated only on ``doFullElectrostatics``, never on ``pressureProfileFreq``, so it evaluates on every full-electrostatics step no matter how rarely profiles are sampled — measured at 14–28× the per-step cost, which would turn example 16's 128000-step sampling stage into 5–11 days per replica.

Replaying a trajectory evaluates once per **sampled frame** instead of once per step, so the cost falls by the sampling stride.  Example 16 samples every 100 steps, so the reconstruction costs about 1% of what computing it inline would have: hours rather than days.

Usage
~~~~~

Hand it the production config the trajectory came from, and the trajectory:

.. code-block:: bash

   $ pestifer pressure-profile-ewald --namd-config prod.namd --dcd prod.dcd --slabs 20 --np 8

The production config is read for its force field, cutoffs, PME settings and restraints, all of which are reused verbatim so that the profile is computed with the same forces that generated the trajectory.  Run control, trajectory output, thermostats and barostats are stripped — nothing is integrated — and the command reports what it removed.

Two configs and two logs are written to ``--workdir``, along with the results:

- ``<prefix>-pressure-profile.csv`` — the frame-averaged profile: slab center :math:`z`, then the three components of the real-space, reciprocal and complete profiles, plus the standard error of the complete one.
- ``<prefix>-pressure-profile.png`` — three panels (lateral, normal, isotropic), each showing both contributions and their sum, so the size of the reciprocal term is visible rather than asserted.

The trajectory must carry unit-cell data — it must have been written with ``DCDUnitCell yes``.  The cell is then read per frame, so a trajectory from an NPT or NPgT run is handled correctly.

Self-check
~~~~~~~~~~

Every run reports how well the reconstruction agrees with NAMD's *own* total pressure.  The real-space pass computes complete forces even though its profile omits the reciprocal term, so the ``PRESSURE`` column of its energy output is an independent reference for the slab-averaged profile:

.. code-block:: text

   INFO> agreement with NAMD's own reported pressure (mean |slab-average - PRESSURE|, over 10 frames):
   INFO>     real space only :     212.73 bar
   INFO>     reconstructed   :       1.05 bar

The real-space half alone is off by hundreds of bar; the reconstruction agrees to about one.  If the reconstructed number is ever the larger of the two, the command warns and the profile should not be trusted.

Options
~~~~~~~

- ``--namd-config`` — the production NAMD config the trajectory came from (required).
- ``--dcd`` — the trajectory to replay (required); must carry unit-cell data.
- ``--slabs`` — ``pressureProfileSlabs`` for both passes (default ``20``).  Both passes necessarily use the same value; summing different slab griddings would be meaningless.
- ``--ewald-grid`` — the ``pressureProfileEwald`` grid, as three integers (default ``10 10 10``, NAMD's own default).  This is a convergence parameter: raise it and re-run to confirm the profile is insensitive to it.
- ``--stride`` — replay every Nth frame (default ``1``).
- ``--temperature`` — temperature used to initialize velocities; affects only the kinetic term (default ``300``).
- ``--workdir`` — directory for the generated configs, logs and output (default ``.``).
- ``--prefix`` — basename for generated files (default ``ppreplay``).
- ``--namd`` — NAMD executable (default ``namd3``).
- ``--np`` — processors for NAMD, passed as ``+pN`` (default ``1``).
- ``--title`` — plot title.
- ``--figsize`` — figure size in inches (default ``15 5``).

A note on frame alignment
~~~~~~~~~~~~~~~~~~~~~~~~~

Each replayed frame is evaluated at timestep 0, so every ``PRESSUREPROFILE:`` record carries the same step number and the frames cannot be matched up by it.  The two passes are aligned **positionally**, and the command refuses to combine passes whose frame or slab counts differ rather than silently pairing the wrong frames.
