.. _subs_mdplot:

mdplot
-------

The ``mdplot`` subcommand generates plots from NAMD log and XST files without requiring a full pestifer configuration file.
It is the standalone equivalent of the ``mdplot`` task used inside a ``pestifer build`` pipeline.

.. code-block:: bash

   $ pestifer mdplot [-h] [--logs LOGS [LOGS ...]] [--basename BASENAME]
            [--figsize FIGSIZE FIGSIZE]
            [--timeseries QUANTITY [QUANTITY ...]]
            [--timecoseries QUANTITY [QUANTITY ...]]
            [--profiles [PROFILE ...]]
            [--profiles-per-block N]
            [--colormap COLORMAP] [--colormap-direction {1,-1}]

Arguments
~~~~~~~~~

``--logs``
   One or more NAMD log files in chronological order.  When multiple logs are supplied they are concatenated before plotting, so runs that were broken into segments (e.g. a restart sequence) are treated as a single continuous trajectory.

``--basename``
   Base name for all output image and CSV files (default: ``mdplot``).

``--figsize``
   Figure width and height in inches (default: ``9 6``).

``--timeseries``
   One or more scalar quantities to plot as individual time-series panels (default: ``density``).  Each quantity produces a separate figure.  Any column appearing in a NAMD ``ENERGY:`` output line or in an XST file is accepted (e.g. ``TOTAL``, ``TEMP``, ``PRESSURE``, ``a_x``).  ``density`` is a pestifer-computed quantity derived from the total system mass (parsed from the log) divided by the ``VOLUME`` column at each step; it is only available for periodic simulations where NAMD reports ``VOLUME``.

``--timecoseries``
   One or more quantities to overlay on a *single* panel.  Use this when you want to compare quantities on the same axes (e.g. ``a_x b_y c_z`` to see all three cell dimensions together).

``--profiles``
   Zero or more profile quantities to plot as a function of position along the z-axis (e.g. ``pressure``).

``--profiles-per-block``
   Number of saved frames to average per profile block (default: ``100``).

``--colormap``
   Matplotlib colormap name used when multiple traces appear on one panel (default: ``viridis``).

``--colormap-direction``
   Colormap direction: ``1`` for normal, ``-1`` for reversed (default: ``1``).

CSV output
~~~~~~~~~~

Providing ``--logs`` triggers two rounds of CSV writing:

1. **Per-log CSVs** — immediately after each log file is parsed, one CSV is written per dataframe type found in that log (``energy``, ``xst``, ``pressureprofile``, …).
2. **Merged CSV** — after all per-log dataframes are concatenated, the combined result is written as ``{basename}-{key}.csv`` (e.g., ``mdplot-energy.csv``).  This is the file that captures the full time series across all supplied log files.

Time-series plot features
~~~~~~~~~~~~~~~~~~~~~~~~~

Each time-series panel uses **simulation time** (in ps or ns, chosen automatically by magnitude) on the bottom x-axis.
The corresponding integer NAMD timestep is shown on the **top spine** as a secondary axis.
Simulation time is computed from the ``TIMESTEP`` and ``TS`` columns in the log, so chained runs with different timestep sizes are handled correctly.

Unit labels for all standard NAMD ``ENERGY:`` output columns (``BOND``, ``ANGLE``, ``DIHED``, ``VDW``, ``ELECT``, ``TOTAL``, ``TEMP``, ``PRESSURE``, ``VOLUME``, etc.) are inferred automatically from NAMD defaults (kcal/mol, K, bar, Å³).

For energetic quantities whose maximum absolute value exceeds 1000 kcal/mol, the axis is automatically scaled by 1/1000 and labeled **1000 kcal/mol** to avoid unwieldy tick magnitudes.

Example
~~~~~~~

Plot total energy and temperature from two chained log files:

.. code-block:: bash

   $ pestifer mdplot --logs run1.log run2.log --timeseries TOTAL TEMP --basename analysis