.. _subs_runtasks_solvate:

solvate 
-------

This task simply performs a standard ``psfgen`` run that invokes the ``solvate`` and ``autoionize`` VMD plugins to make a system solvated by TIP3P water.  It is typically specified using its name only.  A prototypical usage example is below:

.. code-block:: yaml

   tasks:
     - psfgen:
         source:
          id: 6pti
     - md:
         ensemble: minimize
     - md:
         ensemble: NVT
     - solvate:
     - md:
          ensemble: minimize
     - md:
         ensemble: NVT
    
After minimizing and thermalizing (300 K) a newly built system, pestifer will then solvate it, and then it performs a second minimization and thermalization.

