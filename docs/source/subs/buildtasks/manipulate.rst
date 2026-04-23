.. _subs_buildtasks_manipulate:

manipulate
----------

The ``manipulate`` task applies **coordinate-only** modifications to the
pipeline system without changing its topology.  It is intended for rigid-body
operations that reposition atoms — translations, rotations, and alignments —
that need to happen after the PSF has been built (e.g. before solvation or
merging with another system).

Modifications are specified under a ``mods`` key as a dictionary of coormod
types, each containing a list of directives.  Multiple coormod types may be
combined in a single ``manipulate`` task; they are applied in the order they
appear.

Supported coormod types
+++++++++++++++++++++++

.. toctree::
   :maxdepth: 1

   mods/crotations
   mods/transrot
   mods/align

Example
+++++++

The following example orients the system along Z via a backbone alignment to
a reference PDB, then applies an additional translation:

.. code-block:: yaml

   tasks:
     - manipulate:
         mods:
           align:
             - ref_pdb: reference.pdb
               mobile_sel: "backbone"
           transrot:
             - TRANS,0.0,0.0,15.0

For a full description of all parameters see :ref:`config_ref tasks manipulate`.
