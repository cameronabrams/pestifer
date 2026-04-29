.. _installation:

Installation
============

Software Prerequisites
----------------------

The commands ``vmd`` (v. 1.9.4 or better; v. 2.0.0 or better recommended), ``charmrun``, ``catdcd`` (v. 5.2 required), and ``namd3`` (v. 3.0.2 recommended) should be in your path.  By default, Pestifer expects your GPU-enabled ``namd3`` to be in your path as ``namd3gpu``.  To build membrane systems using Pestifer's ``packmol`` integration, you must have ``packmol`` (v. 20.14.3 or better) in your path too.

.. note::

   ``catdcd`` version 5.2 or later is required.  Earlier versions silently drop insertion codes from residue identifiers when reading and writing DCD trajectory files.  Pestifer uses insertion codes to distinguish residues that share the same sequence number (a common occurrence in antibody structures and other proteins with non-standard numbering), so an older ``catdcd`` will corrupt coordinate data for those systems without any warning.

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
   