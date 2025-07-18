.. _config_ref tasks mdplot:

``mdplot``
==========

This task generates plots from the output of MD tasks.  It can plot energy-like quantities vs time, and profile quantities vs distance along the z-axis. The plots are generated using matplotlib, and are saved as image files. This task can be part of a standard build, or it can be used via the command line (``pestifer mdplot``) to generate plots from existing NAMD logs and xst files.


Single-valued parameters:

  * ``postprocessing``: whether to run the mdplot task as part of a standard build or as a post-processing step (default: False)

  * ``profiles``: List of profile quantities to plot vs distance along z-axis

  * ``profiles-per-block``: Number of profiles to average per block (default: 100)

  * ``histograms``: List of histogram quantities to plot vs distance along z-axis

  * ``legend``: include a legend (default: True)

  * ``grid``: include a grid (default: False)

  * ``basename``: Baseame of image/CSV file(s) to generate (default: myplot)



Container-like parameters:

.. toctree::
   :maxdepth: 1

   mdplot/logs
   mdplot/figsize
   mdplot/timeseries
   mdplot/units


.. raw:: html

   <div class="autogen-footer">
     <p>This page was generated by ycleptic v1.7.0 on 2025-07-17.</p>
   </div>