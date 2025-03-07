.. _config_ref namd:

``namd``
========

These are specifications for NAMD, including default input configuration parameters for various classes of simulations.

Single-valued parameters:

  * ``namd-version``: version of NAMD to run (default: 3)

  * ``processor-type``: which processor type to use (cpu or gpu) (default: cpu)



Container-like parameters:

.. toctree::
   :maxdepth: 1

   namd/gpu-resident
   namd/deprecated3
   namd/generic
   namd/vacuum
   namd/solvated
   namd/thermostat
   namd/barostat
   namd/harmonic


