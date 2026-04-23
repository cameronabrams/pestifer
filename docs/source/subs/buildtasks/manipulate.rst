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
   mods/transfer_coords

Example
+++++++

The following example aligns the system to a reference, then injects
refined loop coordinates from a standalone fragment:

.. code-block:: yaml

   tasks:
     - manipulate:
         mods:
           align:
             - ref_pdb: reference.pdb
               mobile_sel: "backbone"
           transfer_coords:
             - donor_pdb: refined_loop.pdb
               donor_sel: "backbone"
               mobile_sel: "chain A and resid 45 to 60 and backbone"
               align_donor_sel: "backbone"
               align_mobile_sel: "chain A and backbone"

For a full description of all parameters see :ref:`config_ref tasks manipulate`.
