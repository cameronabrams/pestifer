.. _installation:

Installation
============

Software Prerequisites
----------------------

Pestifer requires Python 3.12 or newer.  The commands ``vmd`` (v. 1.9.4 or better; v. 2.0.0 or better recommended), ``charmrun``, ``catdcd`` (v. 5.2 required), and ``namd3`` (v. 3.0.2 recommended) should be in your path.  By default, Pestifer expects your GPU-enabled ``namd3`` to be in your path as ``namd3gpu``.  If your GPU-resident NAMD build has a different executable name, override it by setting the ``paths.namd3gpu`` key in your configuration.

.. note::

   ``catdcd`` version 5.2 or later is required.  Earlier versions silently drop insertion codes from residue identifiers when reading and writing DCD trajectory files.  Pestifer uses insertion codes to distinguish residues that share the same sequence number (a common occurrence in antibody structures and other proteins with non-standard numbering), so an older ``catdcd`` will corrupt coordinate data for those systems without any warning.

Installation
------------

To use Pestifer, install it from `PyPI <https://pypi.org/project/pestifer/>`_:

.. code-block:: console

   $ pip install pestifer

Pestifer is under very active development.  To get the latest version, just update:

.. code-block:: console

   $ pip install -U pestifer

.. note::

   The optional ``ligand-paramgen`` extra is needed only by the
   :ref:`make-ligand-mol2 <subs_make_ligand_mol2>` subcommand:

   .. code-block:: console

      $ pip install pestifer[ligand-paramgen]

   It pulls in ``rdkit`` and ``dimorphite_dl``, and the subcommand additionally shells out to
   `Open Babel <https://openbabel.org>`_ (``obabel``), which is a system package rather than a
   Python one -- install it through your OS package manager.

   **This does not give you ligand parameterization.**  ``make-ligand-mol2`` *prepares* a
   CGenFF-ready mol2 file for each unknown ligand; generating the parameters themselves is done by
   the CGenFF program or the `CGenFF web tool <https://cgenff.com/>`_, neither of which ships with
   pestifer.  Pestifer then :ref:`incorporates the resulting stream file <custom_charmm_stream_file>`
   into a build.

   No distributed example requires this extra.  Examples that contain a parameterized ligand --
   9 and 10, whose ``83G`` drug molecule was parameterized with CGenFF beforehand -- carry the
   resulting stream file as a package resource, so they build with nothing extra installed.

If you want a bleeding-edge, potentially unstable version, or you just want the source code, you can clone it from `GitHub <https://github.com/cameronabrams/pestifer>`_:

.. code-block:: console
   
   $ git clone git@github.com:cameronabrams/pestifer.git
   $ cd pestifer
   $ pip install -e .

Of course, you need to pull changes from the repository periodically to keep up with the latest development.  Commits tagged as releases are typically stable.  Consider forking the repository if you want to make changes or contribute to the project.
   