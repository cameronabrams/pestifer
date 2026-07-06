.. _subs_show_resources:

show-resources
--------------

The ``show-resources`` subcommand provides the user a way to view some of the resources packaged along with ``pestifer``.

To see as list of examples provided, use

.. code-block:: console

    $ pestifer show-resources examples

To see a list of all PDB input files suitable for ``make_membrane_system`` tasks, use

.. code-block:: console

    $ pestifer show-resources charmmff --charmmff pdb

To see their full names:

.. code-block:: console

    $ pestifer show-resources charmmff --charmmff pdb --fullnames

Looking up a residue name
=========================

The ``resname`` mode reports, for one or more residue names (case-insensitive), whether pestifer knows the residue: whether it is defined in the CHARMM topology/stream files (and in which file, as a ``RESI`` residue or a ``PRES`` patch), and whether coordinates for it exist in the built-in PDB repository (used by ``make_membrane_system``).

.. code-block:: console

    $ pestifer show-resources resname POPC CHL1 ASPP FOO

produces, for each name, a short block such as

.. code-block:: text

    POPC
      topology: residue (RESI) defined in top_all36_lipid.rtf (segtype: lipid)
      PDB repo: coordinates available (10 conformers) -- "3-palmitoyl-2-oleoyl-D-glycero-1-Phosphatidylcholine"

    ASPP
      topology: patch (PRES) defined in top_all36_prot.rtf
      PDB repo: no coordinates in the built-in PDB repository

    FOO
      topology: not found in any CHARMM topology or stream file
      PDB repo: no coordinates in the built-in PDB repository

This is handy when deciding whether a residue can be built directly (defined in the force field) and, for lipids, whether it can be placed by ``make_membrane_system`` (coordinates in the PDB repository).  For a residue that has coordinates, the block also reports its charge and (for lipids) its head-to-tail length.  Use ``--charmmff-release`` to query a specific CHARMM force-field release and ``--user-pdbcollection`` to include additional PDB collections in the repository search.

To discover residue names rather than look up an exact one, use ``--contains`` to list every known residue name (from both the topologies and the PDB repository) that contains a substring:

.. code-block:: console

    $ pestifer show-resources resname --contains POPC
    2 residue name(s) contain "POPC":
      POPC     RESI lipid     pdb: 10
      POPCE    RESI lipid     pdb: -

Overview
========

Run ``show-resources`` with no argument for a one-line summary of each resource type:

.. code-block:: console

    $ pestifer show-resources
