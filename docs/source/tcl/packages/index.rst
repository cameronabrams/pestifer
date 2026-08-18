.. _tcl-packages:

Tcl packages
============

These packages are used by Pestifer for certain system preparation tasks.  They are not intended
to be run directly by users, but rather are used by Pestifer's Tcl API.  The ``PestiferUtil``
package is loaded automatically when pestifer's Tcl library is initialized (see
:ref:`subs_setup_vmd`); the others must be loaded explicitly with ``package require``.

For example, if you want to use the ``brot`` procedure from the PestiferCRot package, you would
include the following lines at the beginning of your Tcl script:

.. code-block:: tcl

    package require PestiferCRot
    namespace import PestiferCRot::*

Sources for each package are included below for reference.

In use
------

.. tclscript:: pestifer/resources/tcl/pkg/pestifer/util.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/crot.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/environ.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/ionize.tcl

Retained for reference
----------------------

Nothing in pestifer loads the packages below: no Python module and no other Tcl script requires
them.  They ship with pestifer and remain loadable by hand from your own VMD scripts, but they
are not part of any pestifer workflow and are not maintained.  Each carries the same note at the
top of its source.

.. tclscript:: pestifer/resources/tcl/pkg/pestifer/autools.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/axes.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/declash.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/getlinks.tcl
.. tclscript:: pestifer/resources/tcl/pkg/pestifer/multimer.tcl

Third-party packages
--------------------

Pestifer no longer distributes any third-party Tcl package.  The ``La`` (linear algebra) and
``Orient`` packages it used to ship are required only by the retired ``bilayer_orient.tcl``, and
all three now live under ``tests/unit/test_tasks/fixtures/tcl_oracle/`` as the oracle for the
regression test that pins the Python orientation path to the coordinates the Tcl path produced.
The test tree is excluded from both the wheel and the sdist, so none of it reaches an installed
pestifer.
