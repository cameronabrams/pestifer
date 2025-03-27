.. _subs_runtasks_psfgen_mods_transrot:

transrot
--------

A ``transrot`` directive is used to specify global rotations of the entire system (usually before solvation).  There are 2 types of ``transrot`` directives:

1. ``TRANS``: translates the entire system by a specified amount in the x, y, and z directions.  The syntax for specifying a ``TRANS`` directive is ``TRANS,<x>,<y>,<z>``.  ``<x>``, ``<y>``, and ``<z>`` are the amounts of translation in the x, y, and z directions, respectively.  For example, ``TRANS,10.0,0.0,0.0`` would translate the entire system by 10.0 Angstroms in the x direction.
2. ``ROT``: rotates the entire system by a specified amount about the x, y, or z axes.  The syntax for specifying a ``ROT`` directive is ``ROT,<axis>,<angle>``.  ``<axis>`` is the axis about which to rotate (``x``, ``y``, or ``z``), and ``<angle>`` is the angle of rotation in degrees.  For example, ``ROT,x,90`` would rotate the entire system by 90 degrees about the x axis.