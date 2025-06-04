follow-namd-log
---------------

When run interactively, ``pestifer`` automatically shows progress bars for any NAMD run that it starts.  If you want to see the progress of a NAMD run that was started *outside* of pestifer, you can use the ``follow-namd-log`` command.

.. code-block:: console

   $ pestifer follow-namd-log <logname>

Here, ``<logname>`` is the name of the NAMD log file you want to follow.  Pestifer will print a progress bar to the screen as it reads the log file.  The progress bar will update every second or so, and will stop when the NAMD run is finished.  If you want to stop following the log file before it finishes, you can use Ctrl-C to interrupt it.
