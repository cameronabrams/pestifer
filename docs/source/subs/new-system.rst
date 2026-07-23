.. _subs new-system:

new-system
-----------------

Pestifer provides a subcommand ``new-system`` as a convenient way to generate a bare-bones system configuration file for a new simulation.  This is useful when you want to start a new simulation from scratch.

``pestifer new-system`` accepts a single argument, which is intepreted as either the PDB ID or AlphaFold ID of the new system.

For example, to create a new system configuration file for the PDB ID 1abc, you would run:

.. code-block:: bash

   $ pestifer new-system 1abc

This will create a new file named ``1abc.yaml`` in the current directory, containing a basic configuration for the system that, by default, has only a ``fetch`` task and a single ``psfgen`` task:

.. code-block:: yaml

    title: New template pestifer input configuration for PDB ID 1abc
    tasks:
      - fetch:
          sourceID: 1abc
      - psfgen:
           
Using this config is helpful in getting over the first hurdle of a new system, which is making sure pestifer can write a good, albeit minimal, ``psfgen`` script.

The ``new-system`` subcommand bases its output on the example file for example 1, the simple :ref:`BPTI system <example bpti1>` .

There are several options to ``new-system``:

``--full``
+++++++++++

Including this flag option will add solvation, minimization, and equilibration tasks to the configuration file.  For example:

.. code-block:: bash

   $ pestifer new-system 1abc --full

This will create a new file named ``1abc.yaml`` in the current directory, containing a more complete configuration for the system, including tasks for solvation, minimization, and equilibration.

``--inspect``
+++++++++++++

By default ``new-system`` never looks at the actual structure -- it just fills in a template.  The ``--inspect`` flag fetches the structure and reads out the sequence features you would otherwise have to find by hand in the header, then annotates the generated config with the corresponding ready-to-use mod stubs:

.. code-block:: bash

   $ pestifer new-system 4zmj --inspect

It discovers:

- **Missing residues** (``REMARK 465`` / mmCIF ``pdbx_unobs_or_zero_occ_residues``), grouped into contiguous runs and classified as an *interior* gap (both anchors resolved) or an *N-/C-terminal tail* (a free end);
- **Engineered mutations / conflicts** (``SEQADV`` / mmCIF ``struct_ref_seq_dif``) -- residues that differ from the sequence database;
- **Cloning artifacts / expression tags** (``SEQADV`` typed *cloning* / *expression tag*).

Interior missing loops are built and closed automatically, so when any are present ``--inspect`` also adds an active ``ligate`` task to the config.  Every other finding is emitted as a **commented** stub under the ``psfgen`` task -- with the correct nesting shown (a ``sequence:`` block under ``source:``, a ``mods:`` block as a *sibling* of ``source:``) -- so nothing structural is changed until you uncomment the block(s) you want.  For example, the generated ``psfgen`` task for 4zmj includes:

.. code-block:: yaml

   - psfgen:
       source:
         biological_assembly: 1
       # === pestifer new-system --inspect: detected sequence features ===
       # ...
       # To BUILD any as a modeled tail (dropped by default), add under `source:` -
       #   sequence:
       #     terminal_tails:
       #       n: [B, G]
       # To REVERT specific residues, add a `mods:` block (a SIBLING of `source:`) -
       #   mods:
       #     mutations: [B:PRO,559,ILE, G:ASN,332,THR]
       # To EXCISE cloning artifacts / tags, add a `mods:` block -
       #   mods:
       #     deletions: [G:508-513]
   - ligate:
       method: ccd

The inspection is header-based; it does not (yet) align against a canonical UniProt sequence, so a substitution not recorded in ``SEQADV`` will not be flagged.

``--output <filename>``
++++++++++++++++++++++++

This option allows you to specify the name of the output file for the new system configuration.  For example:

.. code-block:: bash

   $ pestifer new-system 1abc --output my_custom_config.yaml

This will create a new file named ``my_custom_config.yaml`` in the current directory, containing the same configuration as ``1abc.yaml``, but with the specified filename.

``--title <title>``
+++++++++++++++++++++++

This option allows you to specify a custom title for the new system configuration.  For example:

.. code-block:: bash

   $ pestifer new-system 1abc --title "My Custom Title"

This will create a new file named ``1abc.yaml`` in the current directory with the specified title.
