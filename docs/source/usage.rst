.. _usage:

Usage
=====

Pestifer's main functionality is available on the command line.  You can also access Pestifer's library of useful TcL/VMD procs from your own VMD scripts.

Command-line Usage
------------------

Installation of Pestifer gives access to the ``pestifer`` command.  The general syntax for invoking ``pestifer`` is

.. code-block:: bash

   $ pestifer <subcommand> <options>

Help with any subcommand can be obtained via

.. code-block:: bash

   $ pestifer <subcommand> --help

General options for all subcommands:

.. code-block:: bash

  --banner, --no-banner             enable or disable the banner
  --log-level {info,debug,warning}  logging level (default: debug)
  --log-file LOG_FILE               log file (default: subcommand-specific)
  --version                         show the pestifer version and exit

Subcommands
-----------

Commands are grouped below by **what you hand them** -- a configuration, a structure, a NAMD run's
output, or nothing at all -- rather than by what they do, so an unfamiliar command can be placed
without reading the whole list.  ``pestifer --help`` prints them under the same headings and in
the same order.

Build a system
~~~~~~~~~~~~~~

The core of pestifer: turn a description of a system into files NAMD can run.

.. toctree::
   :maxdepth: 1

   subs/build
   subs/build-example
   subs/new-system

The ``fetch-example`` subcommand (:ref:`fetch-example <sub_fetch_example>`) is documented alongside
:ref:`build-example <sub_build_example>`; it copies an example's YAML config file to the current
directory without building it.

Configure a build
~~~~~~~~~~~~~~~~~

What you can ask for, and what is available to ask for.  These take no input of yours.

.. toctree::
   :maxdepth: 1

   subs/config-help
   subs/config-default
   subs/show-resources

Work with structures
~~~~~~~~~~~~~~~~~~~~

These take a structure -- a PSF/PDB, a trajectory, a ligand -- rather than a configuration.

.. toctree::
   :maxdepth: 1

   subs/desolvate
   subs/make-ligand-mol2
   subs/make-pdb-collection

After the run
~~~~~~~~~~~~~

These read what NAMD wrote.  All but ``report-methods`` work on **any** NAMD run, including one
pestifer never managed -- they need only the log, the ``.xst`` and the checkpoint files.
``report-methods`` is the exception: it reads the ``run-record.json`` that a pestifer build leaves
behind, so it applies to pestifer's own builds only.

``pressure-profile-ewald`` differs in the other direction: it works on any run, but it is the only
command in this group that *runs NAMD* rather than only reading its output, because the quantity it
reports is one NAMD will not report in a single run.

.. toctree::
   :maxdepth: 1

   subs/mdplot
   subs/density-profile
   subs/pressure-profile-ewald
   subs/follow-namd-log
   subs/make-namd-restart
   subs/report-methods

Manage the installation
~~~~~~~~~~~~~~~~~~~~~~~

These act on pestifer itself rather than on any system of yours.

.. toctree::
   :maxdepth: 1

   subs/cache
   subs/wheretcl
   subs/setup-vmd
   subs/setup-claude
   subs/modify-package

``modify-package`` is a maintainer tool and is registered only when pestifer is run from a source
checkout; it does not appear for a pip-installed pestifer.

.. _custom_charmm_stream_file:

Using a custom CHARMM stream file for a small molecule
------------------------------------------------------

If you have a CHARMM-format stream (``.str``) file for a small molecule -- for
example, a CGenFF stream you generated for a ligand that is not covered by the
standard CHARMM36 release -- pestifer can pull it into a build without any code
or psfgen-script editing.  Add a ``charmmff.user_custom`` block to your build
YAML pointing at the directory that contains the file:

.. code-block:: yaml

   charmmff:
     user_custom:
       searchpath:
         - /path/to/my_ligand_params/

   tasks:
     # ...your normal task list...

At build time, pestifer walks each directory in ``searchpath`` and registers
every CHARMM-format file (``.rtf``, ``.top``, ``.str``, ``.prm``) it finds.
Each topology or stream file is scanned for ``RESI`` and ``PRES`` records, so
the residues and patches it defines become resolvable by name.  When the
``psfgen`` task encounters a residue from your input PDB whose definition lives
in one of these files, the file is automatically copied into the run directory,
emitted as a ``topology <basename>`` directive in the psfgen script, and
registered in the pipeline so its parameter section is also included in the
NAMD parameter file list.

Each residue picked up this way is also classified for the molecule parser.
By default a new residue is treated as a ``ligand`` so it lands in its own
segment.  To override the default, add a ``segtypes`` mapping under
``user_custom``:

.. code-block:: yaml

   charmmff:
     user_custom:
       searchpath:
         - /path/to/my_ligand_params/
       segtypes:
         cofactor: [NADX]

A few things worth knowing:

* A single ``.str`` file containing both ``read rtf … END`` and
  ``read para … END`` blocks is sufficient.  You do not need to also list it
  under ``rtf:`` or ``prm:``.
* The per-extension lists (``rtf:``, ``str:``, ``prm:``) take **basenames only**
  -- any path separator there is rejected.  Directories go under
  ``searchpath:``.
* ``charmmff.user_custom`` is for files you supply.  Do not put them under
  ``charmmff.custom``, which is reserved for files that ship with pestifer.
* Your PDB must use the residue name as it appears in the ``RESI`` card of
  your stream file (e.g. if the stream defines ``RESI LIG``, use ``LIG`` in
  the input PDB).
* Atom names in the input PDB must match the atom names declared inside the
  ``RESI`` block.  RCSB-issued PDBs typically use names that differ from
  CGenFF-generated streams; when they disagree, psfgen builds the topology
  but cannot import coordinates, so each unmatched atom ends up at
  ``(0, 0, 0)``.  Use ``psfgen.aliases.atom`` to bridge the two schemes, e.g.

  .. code-block:: yaml

     psfgen:
       aliases:
         atom:
           - LIG PA P1
           - LIG O1A O1
           - LIG O2A O2

See :ref:`config_ref charmmff user_custom` for the full reference.

.. seealso::

   That ``searchpath`` route keeps the definition local to one build.  To make a
   residue available to *every* build instead, it has to go into pestifer's own
   bundled force field; see :ref:`charmmff_customizations`.

.. _use in vmd scripts:

Use in VMD scripts
------------------

Pestifer provides a library of useful Tcl/VMD procs and VMD macro definitions.
To make them available in every VMD session, run once (from the conda environment
in which pestifer is installed):

.. code-block:: bash

   $ pestifer setup-vmd

This writes ``~/.pestifer/vmd_init.tcl`` with the absolute path to pestifer's Tcl
library hardcoded, and appends a one-line source hook to ``~/.vmdrc``:

.. code-block:: tcl

   # pestifer setup-vmd
   if {[file exists ~/.pestifer/vmd_init.tcl]} { source ~/.pestifer/vmd_init.tcl }

That single line is the only thing that ever needs to be in ``~/.vmdrc`` — it does
not change between pestifer versions.  After upgrading pestifer or switching conda
environments, simply re-run ``pestifer setup-vmd`` to regenerate
``~/.pestifer/vmd_init.tcl`` with the updated path.

.. important::

   Run ``pestifer setup-vmd`` from the shell environment (or conda environment)
   in which pestifer is installed, so that the correct installation path is
   recorded.  VMD itself can be launched from any environment — no dynamic
   ``exec pestifer`` call is made at VMD startup.

Once set up, VMD automatically loads pestifer's Tcl library on startup, including
extended definitions of the ``glycan`` and ``lipid`` atomselect macros and all
PestiferUtil procedures.
