"""
A package for handling building of molecules and generation of their descriptions in psfgen scripts.

The first main task is transfer information from the PDB or mmCIF file.  Since these files by definition define the asymmetric unit of a structure, the :mod:`pestifer.molecule.asymmetricunit` module is responsible for extracting and representing this information in a way that can be used by the rest of the Pestifer framework.  Using the selected biological unit, the :mod:`pestifer.molecule.molecule` module is responsible for building the full structure, including any missing atoms or residues, and generating the necessary files for further processing.

"""