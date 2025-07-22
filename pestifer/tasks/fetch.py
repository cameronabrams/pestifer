# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`InitiateTask` class for initiating a pestifer build.
"""

from ..core.basetask import BaseTask
from ..core.artifacts import PDBFile, CIFFile
from pidibble.pdbparse import PDBParser
import os
import logging
logger=logging.getLogger(__name__)
class FetchTask(BaseTask):
    """
    A class for fetching initial structures from various sources.
    """
    yaml_header='fetch'
    """
    YAML header for the FetchTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """
    _artifact_name = 'base_coordinates'
    def do(self):
        """ Execute the fetch task. """
        self.log_message('initiated')
        source=self.specs.get('source','pdb')
        sourceID=self.specs.get('sourceID','')
        source_format=self.specs.get('source_format','pdb') # desired format for the source file, default is 'pdb'
        if source_format not in ['pdb', 'cif']:
            raise ValueError(f"Unsupported source format: {source_format}. Expected 'pdb' or 'cif'.")
        if source=='pdb':
            if source_format=='pdb':
                # parse the PDB file and register it in the pipeline context
                pdb_file=PDBParser(PDBcode=sourceID).fetch()
                if pdb_file:
                    self.register_current_artifact(PDBFile(sourceID, key=self._artifact_name))
                else:
                    raise ValueError(f"Could not fetch PDB file for sourceID: {sourceID}")
            else:
                cif_file=PDBParser(PDBcode=sourceID, input_format='mmCIF').fetch()
                if cif_file:
                    self.register_current_artifact(CIFFile(sourceID, key=self._artifact_name))
                else:
                    raise ValueError(f"Could not fetch CIF file for sourceID: {sourceID}")
        elif source=='alphafold':
            pdb_file=PDBParser(alphafold=sourceID).fetch()
            if pdb_file:
                self.register_current_artifact(PDBFile(sourceID, key=self._artifact_name))
            else:
                raise ValueError(f"Could not fetch CIF file for sourceID: {sourceID}")
        elif source=='alphafold':
            pdb_file=PDBParser(alphafold=sourceID).fetch()
            if pdb_file:
                self.register_current_artifact(PDBFile(sourceID, key=self._artifact_name))
            else:
                raise ValueError(f"Could not fetch CIF file for sourceID: {sourceID}")
        elif source=='alphafold':
            pdb_file=PDBParser(alphafold=sourceID).fetch()
            if pdb_file:
                self.register_current_artifact(PDBFile(sourceID, key=self._artifact_name))
            else:
                raise ValueError(f"Could not fetch AlphaFold PDB file for sourceID: {sourceID}")
        elif source=='local': # expect <sourceID>.pdb or <sourceID>.cif
            local_pdb=f'{sourceID}.pdb'
            local_cif=f'{sourceID}.cif'
            if source_format=='pdb' and os.path.exists(local_pdb):
                self.register_current_artifact(PDBFile(sourceID, key=self._artifact_name))
            elif os.path.exists(local_cif):
                self.register_current_artifact(CIFFile(sourceID, key=self._artifact_name))
            else:
                raise FileNotFoundError(f"Neither {local_pdb} nor {local_cif} found.")
        else:
            raise ValueError(f"Unsupported source type: {source}")
        self.log_message('complete')
        logger.debug(f'{self.ctx.context_to_string()}')
        return self.result