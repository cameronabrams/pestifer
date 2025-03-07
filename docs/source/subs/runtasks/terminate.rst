.. _subs_runtasks_terminate:

terminate 
---------

A terminate task can be used to build a package that contains everything needed to launch a namd simulation of the newly built system, including copies of all necessary charmmff parameter files.

For example, if your list of tasks ends with the following ``terminate`` task spec:

.. code-block:: yaml

   - terminate:
      basename: my_system
      package:
        ensemble: NPT
        basename: prod_system
        
then pestifer will generate the output files

* ``my_system.psf`` -- the PSF file
* ``my_system.pdb`` -- the PDB file
* ``my_system.coor`` -- the binary NAMD coordinate file equivalent to the PDB file
* ``my_system.xsc``  -- the XSC file that specifies the system box size and shape

It will also generate the tarball ``prod_system.tgz`` that contains

1. The above listed output files;
2. All necessary charmmff parameter files copied from the default toppar directory and gently edited to allow them to be used by NAMD; and
3. An example NAMD config file that you should review and edit.

Transferring the tarball to your production machine is easy, of course.  You may also be building on your production machine (I do), so moving these files and only these files to a clean production directory is a nice thing to do.