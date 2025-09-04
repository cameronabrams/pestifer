# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`InitiateTask` class for initiating a pestifer build.
"""
import logging
import os

from pidibble.pdbparse import PDBParser

from .basetask import BaseTask

from ..core.artifacts import *

logger = logging.getLogger(__name__)

class FetchTask(BaseTask):
    """
    A class for fetching initial structures from various sources.
    """
    _yaml_header = 'fetch'
    """
    YAML header for the FetchTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """
    _artifact_name = 'base_coordinates'
    
    def do(self):
        """ Execute the fetch task. """
        source: str = self.specs.get('source', 'pdb')
        sourceID: str  = self.specs.get('sourceID', '')
        source_format: str = self.specs.get('source_format', 'pdb') # desired format for the source file, default is 'pdb'
        if source_format not in ['pdb', 'cif']:
            raise ValueError(f"Unsupported source format: {source_format}. Expected 'pdb' or 'cif'.")
        if source == 'pdb':
            if source_format == 'pdb':
                # parse the PDB file and register it in the pipeline context
                pdb_file = PDBParser(PDBcode=sourceID).fetch()
                if pdb_file:
                    self.register(sourceID, key=self._artifact_name, artifact_type=PDBFileArtifact)
                else:
                    raise ValueError(f"Could not fetch PDB file for sourceID: {sourceID}")
            else:
                cif_file = PDBParser(PDBcode=sourceID, input_format='mmCIF').fetch()
                if cif_file:
                    self.register(sourceID, key=self._artifact_name, artifact_type=CIFFileArtifact)
                else:
                    raise ValueError(f"Could not fetch CIF file for sourceID: {sourceID}")
        elif source == 'alphafold':
            pdb_file = PDBParser(alphafold=sourceID).fetch()
            if pdb_file:
                self.register(sourceID, key=self._artifact_name, artifact_type=PDBFileArtifact)
                if Path(f'{sourceID}.json').exists():
                    self.register(sourceID, key='json', artifact_type=JSONFileArtifact)
            else:
                raise ValueError(f"Could not fetch CIF file for sourceID: {sourceID}")
        elif source == 'alphafold':
            pdb_file = PDBParser(alphafold=sourceID).fetch()
            if pdb_file:
                self.register(sourceID, key=self._artifact_name, artifact_type=PDBFileArtifact)
            else:
                raise ValueError(f"Could not fetch CIF file for sourceID: {sourceID}")
        elif source == 'alphafold':
            pdb_file = PDBParser(alphafold=sourceID).fetch()
            if pdb_file:
                self.register(sourceID, key=self._artifact_name, artifact_type=PDBFileArtifact)
            else:
                raise ValueError(f"Could not fetch AlphaFold PDB file for sourceID: {sourceID}")
        elif source == 'local':  # expect <sourceID>.pdb or <sourceID>.cif
            local_pdb = f'{sourceID}.pdb'
            local_cif = f'{sourceID}.cif'
            if source_format == 'pdb' and os.path.exists(local_pdb):
                self.register(sourceID, key=self._artifact_name, artifact_type=PDBFileArtifact)
            elif os.path.exists(local_cif):
                self.register(sourceID, key=self._artifact_name, artifact_type=CIFFileArtifact)
            else:
                raise FileNotFoundError(f"Neither {local_pdb} nor {local_cif} found.")
        else:
            raise ValueError(f"Fetch failed: source='{source}' sourceID='{sourceID}'.")
        return self.result