.. _subs_config_help:

config-help
-----------

Because it imports the class :class:`~ycleptic.yclept.Yclept` from `Ycleptic <https://pypi.org/project/ycleptic/>`_, Pestifer has a built-in, interactive, command-line system for help generating YAML-format input configuration files available through the ``config-help`` subcommand.  

.. code-block:: bash

   $ pestifer config-help

       charmmff ->
       psfgen ->
       namd ->
       title
       paths ->
       tasks ->
       .. up
       ! quit
   pestifer-help:


This command ends at a prompt (``pestifer-help:``) that allows you to drill down into the help system.  The help system is organized around the topics that are allowed in a config file, which appear in the list above.  Any item with an arrow after it can be drilled down into.  Double-dot (``..``) takes you up, and bang (``!``) quits.

Any directive can also be reached directly by naming it on the command line, which prints that one entry and leaves you at the prompt:

.. code-block:: bash

   $ pestifer config-help tasks psfgen source exclude

For a worked walkthrough of exploring a config file this way -- and of how the same content reaches the :ref:`config_ref` -- see :ref:`exploring_config_schema`.
