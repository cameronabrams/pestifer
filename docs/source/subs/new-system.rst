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

- **Biological assemblies** (``REMARK 350`` / mmCIF assembly definitions) -- each one's index, how many copies it generates, and which chains it covers -- so you can pick which to build (``biological_assembly: N``; ``0`` = asymmetric unit only);
- **Chain identities** -- every chain with its dominant segment type and size (e.g. ``protein (462 residues)``, ``glycan (NAG, BMA, MAN)``, ``water``, ``ion (ZN)``) -- so you can omit any via ``exclude:`` (for mmCIF the label chain ids are shown, with the author id in parentheses, since those are the ids ``exclude`` and the assembly definitions use);
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

``--interactive``
+++++++++++++++++

``--interactive`` is the guided counterpart to ``--inspect``: it discovers the same features but, instead of writing commented stubs, walks through them one at a time and writes the choices you make as **active** (uncommented) config.

.. code-block:: bash

   $ pestifer new-system 4zmj --interactive
     Build biological assembly 1 (3 copies)? (No = asymmetric unit only) [Y/n] y
     Omit chain A [glycan (NAG, BMA, MAN)]? [y/N] y
     Omit chain G [protein (462 residues)]? [y/N] n
     Build a modeled tail for chain B N-terminus (512-520, 9 res)? [y/N] y
     Add a `ligate` task to close them? [Y/n] y
     Revert B:PRO559 -> db ILE [engineered mutation]? [y/N] y
     ...
     Excise expression tag in chain G (508-513)? [y/N] y

The resulting ``psfgen`` task then carries a real ``source:`` block (chosen ``biological_assembly:`` and an ``exclude:`` list of omitted chains), a ``sequence:`` block (built tails), a ``mods:`` block (chosen mutation reverts and deletions), and an active ``ligate`` task -- ready to run, with no manual editing.  Pressing Enter accepts the shown default (upper-cased in the ``[y/N]`` / ``[Y/n]`` prompt).

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
