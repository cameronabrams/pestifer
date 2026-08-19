.. _subs_buildtasks_mdplot:

mdplot
------

An ``mdplot`` task turns the NAMD logs and ``.xst`` files produced by the ``md`` tasks earlier in
the build into figures.  It is the most commonly used task in pestifer's examples after the ones
that actually build the system — twenty-six of the twenty-seven distributed examples end with one
— because it is how a build shows its work: whether the box density settled, whether the cell
stopped shrinking, how long the run took.

It reads what the build already wrote, so it adds no simulation time.

.. note::

   There is also an :ref:`mdplot subcommand <subs_mdplot>`, which does the same job from the
   command line against logs you already have, without a configuration file.  Use the task to
   have a build produce its own figures; use the subcommand to re-plot afterwards, or to plot a
   production run pestifer did not manage.

A minimal example
=================

.. code-block:: yaml

   - mdplot:
       timeseries:
         - density

Each entry in ``timeseries`` is a quantity to plot against timestep.  A *nested list* puts several
quantities on one set of axes:

.. code-block:: yaml

   - mdplot:
       basename: solvated
       grid: true
       timeseries:
         - - cpu_time          # these two share a plot
           - wall_time
         - density             # this one gets its own
         - - a_x               # and the three cell edges share another
           - b_y
           - c_z

Quantities are whatever the NAMD log and ``.xst`` provide — ``density``, ``temperature``,
``pressure``, the cell edges ``a_x``/``b_y``/``c_z``, ``cpu_time``, ``wall_time``, and the energy
terms NAMD reports.

Membrane systems
================

For a bilayer, two further outputs are usually wanted.  ``profiles`` plots a quantity along the
box normal rather than against time:

.. code-block:: yaml

   - mdplot:
       timeseries:
         - density
       profiles:
         - pressure

which is how the lateral pressure profile figures in :ref:`example mper-tm symmetric bilayer` are
produced.  ``lipids-per-leaflet`` lets the task report area per lipid rather than raw cell area.

Re-plotting a finished build
============================

``reprocess-logs: true`` with an explicit ``logs`` list makes the task read logs that already
exist rather than the ones this build just produced — useful when a task list is being re-run to
regenerate figures only:

.. code-block:: yaml

   - mdplot:
       reprocess-logs: true
       logs:
         - 00-06-000_md-nvt.log
         - 00-07-000_density_equilibrate-NPT.log
       timeseries:
         - density

The logs must be listed in chronological order; the task concatenates them into one series.

Output
======

Figures are written to ``output_dir`` (default ``mdplots/``) with names derived from ``basename``
(default ``myplot``) and the quantities plotted.  ``summary-table: true`` additionally writes a
table of the plotted quantities' converged values.

The full option list — including ``panels``, ``overlay``, ``histograms``, ``block-average``,
``units``, ``colormap`` and the axis-label overrides — is in the
:ref:`configuration reference <config_ref tasks>`.
