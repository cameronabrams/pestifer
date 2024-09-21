alphafold 
---------

If set, its value signifies an accession code used to download a model from the AlphaFold Protein Structure Database. If either ``prebuilt`` or ``id`` are set, this directive is ignored.

Example:

.. code-block:: yaml

  tasks:
    source:
      alphafold: AF-O15552-F1-v4