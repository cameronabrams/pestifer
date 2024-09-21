Configuration File Reference
============================

Pestifer uses YAML-format configuration input for executing a build.  There are six unique top-level directives in a configuration file: ``title``, ``charmmff``, ``namd2``, ``psfgen``, ``paths``, and ``tasks``.  The only one of these directives that is necessary for a build is ``tasks``.  All others have built-in default values that you can override
if you like.  Apart from ``title``, all other directives are detailed below.

.. toctree::

   config_ref/charmmff 
   config_ref/paths
   config_ref/psfgen
   config_ref/namd2
   config_ref/tasks
