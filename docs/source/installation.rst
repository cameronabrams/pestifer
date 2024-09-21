Installation
============


Software Prequisites
--------------------

The commands ``vmd``, ``charmrun``, and ``namd2`` should be in your path.  If you choose to package the results using the ``topogromacs`` plugin, the ``gmx`` executable should also be in your path.  To build membrane systems using ``pestifer``'s ``packmol`` integration, you must have ``packmol`` (v. 20.14.3 or better).

Installation
------------

To use Pestifer, install it from PyPI:

.. code-block:: console

   $ pip install pestifer

Pestifer is under very active development.  To get the latest version from PyPI, just update:

.. code-block:: console

   $ pip install -U pestifer
