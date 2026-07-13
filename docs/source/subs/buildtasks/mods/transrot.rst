.. _subs_buildtasks_psfgen_mods_transrot:

transrot
--------

A ``transrot`` directive specifies a rigid-body transformation of a fully disconnected fragment of the system (usually the whole system, before solvation).  The header ``rottrans`` is accepted as a synonym for ``transrot``.  There are 4 types of ``transrot`` directives:

1. ``TRANS``: translates the fragment by a specified amount in the x, y, and z directions.  The syntax for specifying a ``TRANS`` directive is ``TRANS,<x>,<y>,<z>``.  ``<x>``, ``<y>``, and ``<z>`` are the amounts of translation in the x, y, and z directions, respectively.  For example, ``TRANS,10.0,0.0,0.0`` would translate the fragment by 10.0 Angstroms in the x direction.
2. ``ROT``: rotates the fragment by a specified amount about the x, y, or z axes.  The syntax for specifying a ``ROT`` directive is ``ROT,<axis>,<angle>``.  ``<axis>`` is the axis about which to rotate (``x``, ``y``, or ``z``), and ``<angle>`` is the angle of rotation in degrees.  For example, ``ROT,x,90`` would rotate the fragment by 90 degrees about the x axis.  Rotations are performed about the fragment's own center of mass.
3. ``ALIGN``: rotates the fragment so that a *source* vector is carried onto a *target* vector.  This is given in dictionary form with ``source`` and ``target`` fields (see below).  The rotation is the minimal (roll-free) rotation about the fragment's own center of mass, so no spurious twist about the target axis is introduced.
4. ``AXISANGLE``: rotates the fragment by an angle about the axis normal to a three-atom angle.  Given in dictionary form with ``axis_atoms`` (three atomselections ``[selI, selJ, selK]``) and ``angle``; the pivot is the center of ``selJ`` and the axis is ``(rI - rJ) x (rJ - rK)``.  This is the rigid-body rotation formerly offered as the ``ANGLEIJK`` :ref:`irotation <subs_buildtasks_psfgen_mods_irotations>` (now deprecated).

Vector-based orientation (``ALIGN``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each of ``source`` and ``target`` is either a literal 3-vector ``[x, y, z]`` or a pair of VMD atomselections ``[selA, selB]`` whose mass-weighted centers define the vector (pointing from ``selA`` to ``selB``).  This is convenient for orienting a structure before membrane embedding: define a spanning axis of the protein with two selections and align it onto the membrane normal ``[0, 0, 1]``.

.. code-block:: yaml

    manipulate:
      mods:
        transrot:
          - movetype: ALIGN
            source: ["resid 10 and name CA", "resid 250 and name CA"]  # a spanning axis
            target: [0, 0, 1]                                          # membrane normal (z)

Selecting a fragment
^^^^^^^^^^^^^^^^^^^^^

By default a ``transrot`` acts on the whole system (``sel: all``).  To transform only a part of the system, give the directive in dictionary form with a ``sel`` field holding a VMD atomselection.  A rigid-body transform is only meaningful for a fragment that is *fully disconnected* from the rest of the system, so a disconnection guard is emitted for any selection other than ``all``: if any selected atom is bonded to an atom outside the selection, the build hard-errors rather than silently stretching the crossing bond.

.. code-block:: yaml

    manipulate:
      mods:
        transrot:
          - movetype: ROT
            axis: z
            angle: 90.0
            sel: segid QQQ      # rotate only the fragment QQQ about its own center of mass

Axis-angle rotation (``AXISANGLE``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``AXISANGLE`` rotates the fragment by ``angle`` degrees about the axis normal to the angle defined by three atoms.  Give it ``axis_atoms``, a list of three VMD atomselections ``[selI, selJ, selK]`` (each normally selecting a single atom); the pivot is the center of ``selJ`` and the rotation axis is ``(rI - rJ) x (rJ - rK)``.  As with the other movetypes, ``sel`` names the fragment to rotate (default ``all``, disconnection-guarded).

.. code-block:: yaml

    manipulate:
      mods:
        transrot:
          - movetype: AXISANGLE
            axis_atoms:
              - "segid A and resid 12 and name N"
              - "segid B and resid 15 and name C"
              - "segid B and resid 18 and name O"
            angle: 180.0
            sel: segid B

This is the rigid-body rotation formerly offered as the ``ANGLEIJK`` :ref:`irotation <subs_buildtasks_psfgen_mods_irotations>`.  ``ANGLEIJK`` still works there for backward compatibility but is deprecated in favor of ``AXISANGLE``, since it is a rigid-body operation rather than an internal-coordinate one.
