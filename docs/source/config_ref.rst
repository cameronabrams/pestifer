.. _config_ref:

Configuration File Reference
============================

Pestifer uses YAML-format configuration input for executing a build.  There are six unique top-level directives in a configuration file:

Single-valued parameter:

  * ``title``: Meaningful title (optional) (default: Pestifer)

Subdirectives:

.. toctree::
   :maxdepth: 1

   config_ref/charmmff
   config_ref/psfgen
   config_ref/namd
   config_ref/paths
   config_ref/tasks


This documentation was created automatically by ``yclept make-doc`` from the Ycleptic package (`pypi <https://pypi.org/project/ycleptic/>`_) using ``pestifer/resources/ycleptic/base.yaml`` as input.