.. _tcl-scripts:

Tcl scripts
===========

These scripts are used by Pestifer for certain system preparation tasks.  They are not intended
to be run directly by users, but rather are used by Pestifer's Tcl API.

In use
------

These two are sourced by the :ref:`make_membrane_system <subs_buildtasks_make_membrane_system>`
task as it builds and embeds a bilayer.

.. tclscript:: pestifer/resources/tcl/scripts/bilayer_embed.tcl
.. tclscript:: pestifer/resources/tcl/scripts/bilayer_patch.tcl

``bilayer_orient.tcl`` was formerly sourced here as well.  Orientation now happens in Python
(:class:`~pestifer.objs.rottrans.RotTrans`, via ``MakeMembraneSystemTask._orientation_align``),
so no pestifer workflow runs it.  It has moved out of the package to
``tests/unit/test_tasks/fixtures/tcl_oracle/``, where it serves as the oracle for the regression
test that pins the Python path to the coordinates the Tcl path produced.

Retained for reference
----------------------

The domain-swap task was retired and is no longer registered as a user-invocable task, so
nothing in a pestifer workflow runs this script.  It is kept so the method is not lost.

.. _tcl-domainswap:
.. tclscript:: pestifer/resources/tcl/scripts/domainswap.tcl
