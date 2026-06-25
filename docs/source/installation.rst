.. _installation:

Installation
============

Software Prerequisites
----------------------

The commands ``vmd`` (v. 1.9.4 or better; v. 2.0.0 or better recommended), ``charmrun``, ``catdcd`` (v. 5.2 required), and ``namd3`` (v. 3.0.2 recommended) should be in your path.  By default, Pestifer expects your GPU-enabled ``namd3`` to be in your path as ``namd3gpu``.  To build membrane systems with the ``packmol`` packer (``bilayer.packer: packmol``, the default), you must have ``packmol`` (v. 20.15.1 or better) in your path too; the ``grid`` packer places lipids directly and needs no ``packmol``.

.. note::

   ``catdcd`` version 5.2 or later is required.  Earlier versions silently drop insertion codes from residue identifiers when reading and writing DCD trajectory files.  Pestifer uses insertion codes to distinguish residues that share the same sequence number (a common occurrence in antibody structures and other proteins with non-standard numbering), so an older ``catdcd`` will corrupt coordinate data for those systems without any warning.

.. note::

   **packmol and conda environments.**  AmberTools ships a ``packmol`` binary as part of ``packmol-memgen``, and any active conda environment that has AmberTools installed will put that (older) ``packmol`` at the front of your ``PATH``.  Pestifer will refuse to run with it because it predates the ``pbc`` keyword.

   Bash aliases (``alias packmol=/usr/local/bin/packmol``) do **not** help here, because pestifer launches ``packmol`` from a non-interactive subprocess where shell aliases are not expanded.

   Recommended best practice:

   1. Download packmol >= 20.15.1 from `the m3g/packmol GitHub repository <https://github.com/m3g/packmol>`_ and build it.
   2. Install the resulting binary somewhere stable outside any conda environment, e.g. ``/usr/local/bin/packmol``.
   3. Tell pestifer to use that binary by setting ``paths.packmol`` in your run YAML:

      .. code-block:: yaml

         paths:
           packmol: /usr/local/bin/packmol

   The ``paths.packmol`` setting is an absolute path that bypasses ``PATH`` resolution entirely, so it cannot be shadowed by conda activation.

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
   