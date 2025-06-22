.. _subs_mdplot:

mdplot
-------

The ``mdplot`` subcommand is used to plot data from molecular dynamics simulations. It can generate various types of plots, such as energy, temperature, and pressure over time, or radial distribution functions.
It is particularly useful for analyzing the results of molecular dynamics simulations and visualizing the behavior of the system over time.

.. code-block:: console

   $ pestifer mdplot [-h] [--logs LOGS [LOGS ...]] [--xsts XSTS [XSTS ...]] 
            [--basename BASENAME] [--savedata SAVEDATA] [--figsize FIGSIZE FIGSIZE]
            [--traces TRACES [TRACES ...]] [--profiles PROFILES [PROFILES ...]]

The ``mdplot`` subcommand mimics the ``mdplot`` task in the ``pestifer run`` command, but it does not require a configuration file. Instead, it allows you to specify the NAMD log and xst files directly on the command line.