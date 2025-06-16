.. _subs_show_resources:

show-resources
--------------

The ``show-resources`` subcommand provides the user a way to view some of the resources packaged along with ``pestifer``.

To see as list of examples provided, use

.. code-block:: console

    $ pestifer show-resources --examples

To see a list of all PDB input files suitable for ``make_membrane_system`` tasks, use

.. code-block:: console

    $ pestifer show-resources --charmmff pdb

To see their full names:

.. code-block:: console

    $ pestifer show-resources --charmmff pdb --fullnames

