.. _installation:

Installation
============


Software Prequisites
--------------------

The commands ``vmd``, ``charmrun``, and ``namd3`` should be in your path.  By default, ``pestifer`` expects your GPU-enabled ``namd3`` to be in your path as ``namd3gpu``.  To build membrane systems using ``pestifer``'s ``packmol`` integration, you must have ``packmol`` (v. 20.14.3 or better).

Installation
------------

To use Pestifer, install it from `PyPI <https://https://pypi.org/project/pestifer/>`_:

.. code-block:: console

   $ pip install pestifer

Pestifer is under very active development.  To get the latest version, just update:

.. code-block:: console

   $ pip install -U pestifer

If you want a bleeding-edge, potentially unstable version, or you just want the source code, you can clone it from `GitHub <https://github.com/cameronabrams/pestifer>`_:

.. code-block:: console
   
   $ git clone git@github.com:cameronabrams/pestifer.git
   $ cd pestifer
   $ pip install -e .

   