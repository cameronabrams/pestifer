.. _installation:

Installation
============

Software Prequisites
--------------------

The commands ``vmd`` (v. 1.9.4 or better), ``charmrun``, ``catdcd``, and ``namd3`` (v. 3.0.2 recommended) should be in your path.  By default, Pestifer expects your GPU-enabled ``namd3`` to be in your path as ``namd3gpu``.  To build membrane systems using Pestifer's ``packmol`` integration, you must have ``packmol`` (v. 20.14.3 or better) in your path too.

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

Of course, you need to pull changes from the repository periodically to keep up with the latest development.  Commits tagged as releases are typically stable.  Consider forking the repository if you want to make changes or contribute to the project.
   