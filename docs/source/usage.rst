Usage
=====

Installation of the ``pestifer`` package gives access to the ``pestifer`` command.  Pestifer expects an input configuration file for any run.

.. code-block:: console

   $ pestifer my_config.yaml

So to use pestifer, one really needs to know how to set up the configuration file.

An example YAML configuration file to build the 6PTI system in the psfgen user manual is below:

.. code-block::yaml

   title: BPTI
   tasks:
   - psfgen:
      source:
         rcsb: 6pti
         biological_assembly: 1
      minimize:
         nminsteps: 1000
      cleanup: True
   - solvate:
      minimize:
        nminsteps: 1000
   - relax:
      temperature: 300
      pressure: 1
      nsteps: 1000
      dcdfreq: 100   
         
Pestifer is instructed to execute a series of ``tasks``.  The first is a ``psfgen`` task, where the experimental coordinates are used with psfgen to generate an initial topology/coordinate file couple. The ``source`` block declares we will use the 6pti entry of the RCSB, and we will generate biological assembly 1 (which is the only one for this file).  We then minimize the system using NAMD2, followed by cleaning up any intermediate files generated duing psfgen.  The second task is solvation and ionization, which is also followed by minimization using NAMD2.  Finally, the system is relaxed at a temperature of 300 K and a pressure of 1 bar.