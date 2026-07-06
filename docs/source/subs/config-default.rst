.. _subs_config_default:

config-default
--------------

The ``config-default`` subcommand prints the default values of pestifer's configuration options.  With no arguments it prints every default; given one or more directive names it restricts the output to those directives.  It is the non-interactive companion to :ref:`config-help <subs_config_help>`, which lets you explore the same configuration schema interactively.

Print all defaults:

.. code-block:: bash

   $ pestifer config-default

Restrict to particular directives (e.g. the defaults for the ``md`` task):

.. code-block:: bash

   $ pestifer config-default tasks md
