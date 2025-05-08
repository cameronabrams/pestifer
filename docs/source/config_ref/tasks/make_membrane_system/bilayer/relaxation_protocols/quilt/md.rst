.. _config_ref tasks make_membrane_system bilayer relaxation_protocols quilt md:

``md``
======

md task for bilayer relaxation

Single-valued parameters:

  * ``ensemble``: Name of ensemble to run in; note that NPT can only be performed on a solvated system (default: minimize)

  * ``nsteps``: Number of MD time steps (or minimization steps if 'ensemble' is 'minimize') (default: 1000)

  * ``temperature``: Temperature (K) for NVT or NPT (default: 300)

  * ``pressure``: Pressure (bar) for NPT; note that NPT can only be performed on a solvated system (default: 1)



Container-like parameter:

.. toctree::
   :maxdepth: 1

   md/addl_parmfiles


Subdirective:

.. toctree::
   :maxdepth: 1

   md/other_parameters


