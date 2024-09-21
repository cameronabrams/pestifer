id 
--

This single-valued directive allows you to specify a 4-character PDB ID of the source, which pestifer can download, or the basename (without the .pdb extension) of a pdb file in the current working directory.

Example:

.. code-block:: yaml

  tasks:
    psfgen:
      id:  1gc1

