.. _config_ref tasks terminate package:

``package``
===========

Task that packages all CHARMM parameter files and pdb/psf/coor/vel/xsc files for a production run

Single-valued parameters:

  * ``basename``: Desired basename for final tarball (default: my_system)

  * ``ensemble``: Name of ensemble to run in; note that NPT can only be peformed on a solvated system (default: NVT)

  * ``nminsteps``: Number of minimization steps (default: 0)

  * ``nsteps``: Number of MD time steps (default: 1000000)

  * ``dcdfreq``: number of time steps between configuration output to DCD file (default: 10000)

  * ``xstfreq``: number of time steps between cell size output to DCD file (default: 10000)

  * ``restartfreq``: number of timesteps between restart checkpoints (default: 10000)

  * ``temperature``: Temperature (K) for NVT or NPT (default: 300)

  * ``pressure``: Pressure (bar) for NPT; note that NPT can only be peformed on a solvated system (default: 1)

  * ``topogromacs``: Set to True if you want to use topogromacs to generate Gromacs-compatible PDB and TOP files (default: False)



Subdirective:

.. toctree::
   :maxdepth: 1

   package/constraints


