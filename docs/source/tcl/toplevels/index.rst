.. _tcl-toplevels:

Top-level scripts
=================

These are Tcl scripts that are always run when VMD is started within a Pestifer workflow.  They can also be run manually from the VMD console using the ``pestifer_init`` Tcl command, provided you have defined that ``proc`` in your own VMD startup script (e.g., ``~/.vmdrc``).  See the :ref:`use in vmd scripts` section of the :ref:`usage` documentation for more information.


.. tclscript:: pestifer/resources/tcl/vmdrc.tcl

.. tclscript:: pestifer/resources/tcl/macros.tcl