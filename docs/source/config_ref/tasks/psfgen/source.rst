``source``
==========

Specifies the source of the initial coordinate file

Single-valued parameters:

  * ``id``: The 4-character PDB ID of the source or the basename of a local coordinate file (PDB or mmCIF format); pestifer will download from the RCSB if a file is not found

  * ``alphafold``: If set, its value signifies an accession code used to download a model from the AlphaFold database. If either prebuilt or id are set, this directive is ignored

  * ``biological_assembly``: integer index of the biological assembly to construct; default is 0, signifying that the asymmetric unit is to be used (default: 0)

  * ``transform_reserves``: dictionary keyed by chain IDs in the asymmetric unit with values that are lists of chain ID's for symmetry-related chains in other units of the biological assembly

  * ``remap_chainIDs``: dictionary mapping of chainIDs in structure input to new chainIDs

  * ``reserialize``: If true, resets the serial numbers of all atoms based on their order in the input structure file -- use with caution!! (default: False)

  * ``model``: model number to use if there are multiple models in the file

  * ``file_format``: either PDB or mmCIF; some entries do not have a PDB-format file.  The main advantage of PDB is that it uses the author-designations for chains by default.  mmCIF is the new "default" format of the PDB. (default: PDB)

  * ``cif_residue_map_file``: name of output file to contain a mapping of CIF chain-residue number to author chain-residue number-insertion



Subdirectives:

.. toctree::
   :maxdepth: 1

   source/prebuilt
   source/exclude
   source/sequence


