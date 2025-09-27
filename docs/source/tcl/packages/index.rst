.. _tcl-packages:

Tcl packages
============

These packages are used by Pestifer for certain system preparation tasks.  They are not intended to be run directly by users, but rather are used by Pestifer's Tcl API.  The procedures they define are available in the Pestifer Tcl API, and can be used in your own VMD scripts.  The PestiferUtil package is automatically loaded by the ``pestifer_init`` procedure.  Others must be loaded explicitly using the ``package require`` command.

For example, if you want to use the ``brot`` procedure from the PestiferCRot package, you would include the following lines at the beginning of your Tcl script:

.. code-block:: tcl

    package require PestiferCRot
    namespace import PestiferCRot::*

Sources for each package are included below for reference.

.. tclscript:: pestifer/resources/tcl/pkg/pestifer/autools.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/axes.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/crot.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/declash.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/environ.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/getlinks.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/multimer.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/util.tcl
