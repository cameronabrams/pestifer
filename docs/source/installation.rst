Installation
============


Software Prequisites
--------------------

The commands ``vmd``, ``charmrun``, and ``namd2`` should be in your path.  If you choose to package the results using the ``topogromacs`` plugin, the ``gmx`` executable should also be in your path.  To build membrane systems using ``pestifer``'s ``packmol_memgen`` integration, you must have the ``ambertools`` package (v. 23.6 or better) installed (preferable via ``conda``).

Installation
------------

To use Pestifer, install it from PyPI:

.. code-block:: console

   $ pip install pestifer

Pestifer is under very active development.  To get the latest version from PyPI, just update:

.. code-block:: console

   $ pip install -U pestifer
