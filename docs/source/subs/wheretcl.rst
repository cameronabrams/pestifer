.. _subs wheretcl:

wheretcl
--------

Pestifer has a pretty handy library of TcL packages.  If you want to peruse the source, pestifer will tell you where to find them:

.. code-block:: console

   $ pestifer wheretcl --pkg-dir

This subcommand is mainly exposed to enable VMD sessions launched *outside* of pestifer to use pestifer's Tcl library; see :ref:`use in vmd scripts`. 