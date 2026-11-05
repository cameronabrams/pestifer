.. _subs_pressure_profile_ewald:

pressure-profile-ewald
----------------------

The ``pressure-profile-ewald`` subcommand completes a NAMD pressure profile.  It replays a production trajectory for the PME reciprocal-space term, which NAMD leaves out of the ordinary profile, and adds it frame by frame to the profile the production run already wrote.

Why the reciprocal term has to be added
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A run with ``pressureProfile on`` reports, in NAMD's own words, the "kinetic, bonded and nonbonded (but not reciprocal space) contributions".  The missing reciprocal term is not a small correction.  On a 200-lipid DMPC patch it is about :math:`-180` bar of the trace and :math:`+33` mN/m of the integrated tension, and it is anisotropic, so it changes the surface tension itself, not only the pressure.

The obvious fix — switching on ``pressureProfileEwald`` as well — does not work, because NAMD makes the two mutually exclusive.  From ``SimParameters.C`` (3.0.3, line 6699):

.. code-block:: c

   // if Ewald is on, only calculate Ewald
   if (pressureProfileEwaldOn)
     pressureProfileOn = 0;

So a run with ``pressureProfileEwald on`` reports the reciprocal contribution **alone**.

Nor is it practical to compute the reciprocal term during the production run.  NAMD registers ``ComputeEwald`` as an ordinary compute gated only on ``doFullElectrostatics``, never on ``pressureProfileFreq``, so it evaluates on every full-electrostatics step no matter how rarely profiles are sampled — measured at 14–28× the per-step cost.  Replaying a trajectory evaluates it once per **sampled frame** instead, so the cost falls by the sampling stride.

What the production run must do
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- ``pressureProfile on``.  Its profile is the real-space half of the result.
- ``fullElectFrequency 1``.  Under multiple time stepping NAMD's ordinary profile already carries part of the reciprocal term, so adding the replayed term counts that part twice.  On the DMPC patch, ``fullElectFrequency 2`` put the sum 70–95 bar per component away from NAMD's own pressure tensor, about +14.5 mN/m of total tension that does not exist.  At ``fullElectFrequency 1`` it agreed to about 1 bar.  No rescaling of the reciprocal term repairs it exactly, so the command refuses an MTS production config unless given ``--accept-mts``.  Pestifer writes ``fullElectFrequency 1`` into any stage it generates with a pressure profile on.
- ``dcdfreq`` a multiple of ``pressureProfileFreq``, with ``DCDUnitCell yes``.  Each DCD frame is matched to the production profile record **at the same MD step**.  The cell is read per frame, so NPT and NPγT trajectories are handled correctly.
- ``outputPressure`` equal to ``pressureProfileFreq``, so that the log's ``PRESSAVG`` tensor averages the same window as each profile record.  This is what the result is checked against.

Usage
~~~~~

.. code-block:: bash

   $ pestifer pressure-profile-ewald --namd-config prod.namd --log prod.log --dcd prod.dcd --np 8

The production config is read for its force field, cutoffs, PME settings and restraints, all of which are reused verbatim so that the reciprocal term is computed with the same forces that generated the trajectory.  Run control, trajectory output, thermostats and barostats are stripped — nothing is integrated — and the command reports what it removed.

The replay config and log are written to ``--workdir``, along with the results:

- ``<prefix>-pressure-profile.csv`` — the frame-averaged profile: slab center :math:`z` (measured from the bottom of the cell, ``origin.z - c/2``), then the three components of the real-space, reciprocal and complete profiles, plus the standard error of the complete one.
- ``<prefix>-pressure-profile.png`` — three panels (lateral, normal, isotropic), each showing both contributions and their sum.

Self-check
~~~~~~~~~~

Every run compares the slab average of the complete profile with the production run's own ``PRESSAVG`` tensor, **component by component**:

.. code-block:: text

   INFO> agreement with the production run's own pressure tensor (frame mean of slab average - PRESSAVG, xx / yy / zz, over 200 frames):
   INFO>     real space only :  +194.14 /  +196.22 /  +152.46 bar
   INFO>     complete        :    -1.06 /    -1.09 /    -0.90 bar

That is the DMPC patch: the real-space half alone is off by about 150–200 bar, and the complete profile agrees to about one.  The command warns, and the profile should not be used, when any component misses by more than 10 bar (or three standard errors of the per-frame difference, if that is larger).  The likely causes, most common first:

- the Ewald grid is too coarse for the box — re-run with a finer ``--ewald-grid`` and see whether the gap shrinks;
- the production config no longer matches the run, or the log is from a different run;
- a NAMD build that writes wrong pressure profiles while its energies are right.  That one misses by thousands of bar, not tens; see below.

The check is per component because a surface tension is the *difference* between the normal and lateral components.  An earlier version of this command checked only the trace, against a pressure its own replay computed, and could see neither an anisotropic error nor one the replay shared.

If the production log has no ``PRESSAVG`` lines, the command says the result was **not** checked rather than passing it.

The Ewald grid
~~~~~~~~~~~~~~

By default the ``pressureProfileEwald`` grid is sized from the first frame's cell at about 3 Å per point.  NAMD's own default, 10 points per axis, is far too coarse for a membrane box.  On example 17's ~115 × 116 × 81 Å cell the reconstruction missed NAMD's pressure by 182 bar at 10³, 53 bar at 20³, 1.2 bar at 30³ and 0.9 bar at 40³.  The cost of the replay grows with the number of grid points, so do not go far past convergence.

Checking a profile from any run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A NAMD build can compute energies, forces, ``PRESSURE`` and the barostat correctly while writing a wrong per-slab profile.  Such a build's trajectories are fine, and nothing in its log looks wrong.  To check a binary you rely on for profiles:

1. Average :math:`(P_{xx} + P_{yy} + P_{zz})/3` over all slabs and all frames.
2. Compare the result with the mean of the ``PRESSURE`` column over the same frames.

For an ordinary ``pressureProfile on`` run, the two should agree to within a few hundred bar; the gap is the reciprocal-space term that profile leaves out.  A gap of thousands of bar, the same on every frame, means the binary is at fault, not the system.  One such build gave a slab average near :math:`-8000` bar against a reported pressure of a few hundred.

Options
~~~~~~~

- ``--namd-config`` — the production NAMD config the trajectory came from (required).
- ``--log`` — the production NAMD log (required): its ``PRESSUREPROFILE`` records and ``PRESSAVG`` tensor.
- ``--dcd`` — the production trajectory (required); must carry unit-cell data.
- ``--ewald-grid`` — the ``pressureProfileEwald`` grid, as three integers (default: sized from the cell at about 3 Å).  Raise it and re-run to confirm the profile is insensitive to it.
- ``--accept-mts`` — proceed even though the production run used ``fullElectFrequency`` above 1.  The result then double-counts part of the reciprocal term.
- ``--slabs`` — optional check of the slab count; the replay always uses the production profile's own.
- ``--stride`` — replay every Nth frame (default ``1``).
- ``--temperature`` — the temperature NAMD needs to start; it does not enter the reciprocal term (default ``300``).
- ``--workdir`` — directory for the generated config, log and output (default ``.``).
- ``--prefix`` — basename for generated files (default ``ppreplay``).
- ``--namd`` — NAMD executable (default ``namd3``).
- ``--np`` — processors for NAMD, passed as ``+pN`` (default ``1``).
- ``--title`` — plot title.
- ``--figsize`` — figure size in inches (default ``15 5``).

Why there is no real-space replay
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Versions before 3.31.0 also replayed the real-space profile and summed the two replays.  That half was wrong.  A replay evaluates each frame at ``run 0`` with one velocity set drawn from ``--temperature`` and reused for every frame, and a single force evaluation carries no rigid-bond constraint virial.  On the DMPC patch the replayed real-space profile differed from the production run's own by tens of bar per component with rigid bonds and by about 485 bar of the trace without them.  Its self-check passed regardless, because it compared the replay with the replay's own ``PRESSURE``, computed from the same velocities.
