.. _subs_setup_vmd:

setup-vmd
---------

The ``setup-vmd`` subcommand configures VMD to load pestifer's Tcl library automatically on startup.

.. code-block:: console

   $ pestifer setup-vmd

This writes ``~/.pestifer/vmd_init.tcl`` with the absolute path to pestifer's Tcl root hardcoded, and appends a one-line source hook to ``~/.vmdrc``:

.. code-block:: tcl

   # pestifer setup-vmd
   if {[file exists ~/.pestifer/vmd_init.tcl]} { source ~/.pestifer/vmd_init.tcl }

The hook in ``~/.vmdrc`` never needs to change between pestifer versions.  After upgrading pestifer or switching conda environments, simply re-run ``pestifer setup-vmd`` to regenerate ``~/.pestifer/vmd_init.tcl`` with the updated path.

.. important::

   Run ``pestifer setup-vmd`` from the conda environment in which pestifer is installed so that the correct installation path is recorded.  VMD itself can be launched from any environment.

Once set up, VMD automatically loads pestifer's Tcl library on startup, including extended definitions of the ``glycan`` and ``lipid`` atomselect macros and all PestiferUtil procedures.
