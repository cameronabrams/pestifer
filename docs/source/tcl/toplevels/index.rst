.. _tcl-toplevels:

Top-level scripts
=================

These are the Tcl scripts that run when VMD is started within a Pestifer workflow.
``vmdrc.tcl`` puts pestifer's packages on VMD's ``auto_path``, loads ``PestiferUtil``, and
sources ``macros.tcl`` to extend VMD's atomselect macros.

To get the same environment in your own interactive VMD sessions, run ``pestifer setup-vmd``
once; see :ref:`subs_setup_vmd` and :ref:`use in vmd scripts`.

.. tclscript:: pestifer/resources/tcl/vmdrc.tcl

.. tclscript:: pestifer/resources/tcl/macros.tcl
