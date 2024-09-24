``md``
======

Parameters controlling a NAMD run

Single-valued parameters:

  * ``vacuum``: if true, namd will run on a single node even if multiple nodes are allocated (default: False)

  * ``ensemble``: Name of ensemble to run in; note that NPT can only be peformed on a solvated system (default: NVT)

  * ``minimize``: Number of minimization steps; ignored if ensemble is not 'minimize' (used with the 'minimize' command) (default: 1000)

  * ``nsteps``: Number of MD time steps; ingored if ensemble is 'minimize' (used with the 'run' command) (default: 2000)

  * ``dcdfreq``: number of time steps between configuration output to DCD file (default: 100)

  * ``xstfreq``: number of time steps between cell size output to DCD file (default: 100)

  * ``temperature``: Temperature (K) for NVT or NPT (default: 300)

  * ``pressure``: Pressure (bar) for NPT; note that NPT can only be peformed on a solvated system (default: 1)

  * ``other_parameters``: key:value pairs for other namd2 configuration file statements



Subdirective:

.. toctree::
   :maxdepth: 1

   md/constraints


