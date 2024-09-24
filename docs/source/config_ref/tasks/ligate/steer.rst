``steer``
=========

Specifies parameters for steered MD to bring C-termini close to their connection points

Single-valued parameters:

  * ``ensemble``: Name of ensemble to run in; note that NPT can only be peformed on a solvated system (default: NVT)

  * ``force_constant``: force constant (kcal/mol/Angstrom^2) (default: 20)

  * ``target_distance``: distance (Angstrom) between carbonyl carbon of loop C-terminus and amide nitrogen of connecting residue that steered MD drives to (default: 2)

  * ``nsteps``: number of MD time steps for steeering (default: 1000)

  * ``dcdfreq``: number of time steps between configuration output to DCD file (default: 100)

  * ``xstfreq``: number of time steps between cell size output to DCD file (only used if ensemble is NPT) (default: 100)

  * ``temperature``: Temperature used in steered MD run (default: 300)

  * ``pressure``: Pressure used in steered MD run (only used if ensemble is NPT) (default: 1)

  * ``receiver_flexible_zone_radius``: size of zone around receiver N-terminal N in which CA's are not fixed during steered MD (default: 0.0)



Subdirective:

.. toctree::
   :maxdepth: 1

   steer/constraints


