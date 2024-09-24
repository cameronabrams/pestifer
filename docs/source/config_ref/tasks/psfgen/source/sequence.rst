``sequence``
============

Parameters controlling sequence modifications based on specifications in the PDB source

Single-valued parameters:

  * ``include_terminal_loops``: Specifies whether or not N- and C- terminal loops (if unresolved) are built in (default: False)

  * ``fix_engineered_mutations``: Specifies whether or not sequence differences w.r.t. a sequence database and which are labeled as "engineered mutation"s are reverted to their database values (default: True)

  * ``fix_conflicts``: Specifies whether or not sequence differences w.r.t. a sequence database and which are labeled as "conflicts"s are reverted to their database values (default: True)



Container-like parameters:

.. toctree::
   :maxdepth: 1

   sequence/build_zero_occupancy_C_termini
   sequence/build_zero_occupancy_N_termini


Subdirectives:

.. toctree::
   :maxdepth: 1

   sequence/loops
   sequence/glycans


