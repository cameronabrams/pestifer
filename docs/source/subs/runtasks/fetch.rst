.. _subs_runtasks_fetch:

fetch
-----

A ``fetch`` task allows you to retrieve a structure file from a remote source, such as the PDB or AlphaFold databases.
This task is typically the first step in a pestifer run, injecting the first representation of the system state (the downloaded structure file) into the workflow.

The fetch task currently supports two sources:

1. **PDB**: The Protein Data Bank (PDB) is a repository for 3D structural data of biological macromolecules. You can specify a PDB ID to fetch a structure from the PDB.

2. **AlphaFold**: AlphaFold is a deep learning model for predicting protein structures. You can specify a UniProt ID to fetch a structure from the AlphaFold database.

A ``fetch`` task has three attributes that can be set in the config file:

1. ``source``: either ``pdb`` or ``alphafold`` (case-insensitive); ``pdb`` by default.
2. ``sourceID``: the ID of the structure in the source database (a PDB ID or a UniProt ID); **this attribute is required**.
3. ``source_format``: either ``pdb`` or ``cif``.  Some structures in the PDB are only available in the new mmCIF/PDBx format; in this case you should specify ``cif`` here.  Default is ``pdb``.

For example, to begin a run using the bovine pancreatic trypsin inhibitor (BPTI) `structure PDB ID 6pti <https://www.rcsb.org/structure/6PTI>`_ from the PDB, you would specify the following in your config file:

.. code-block:: yaml
        
    tasks:
    - fetch:
        sourceID: 6pti

If the structure file has already been downloaded (i.e., Pestifer detects it in your CWD), no download is performed.