# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`InitiateTask` class for initiating a pestifer build.
"""

from ..core.basetask import BaseTask
from ..core.pipeline import PDBFile, CIFFile
from pidibble.pdbparse import PDBParser
import os

class FetchTask(BaseTask):
    """
    FetchTask class for fetching data in a pestifer build.
    This class inherits from the :class:`BaseTask` class and is used to set up the initial state of the build.
    It prepares the environment, initializes necessary configurations, and sets up the context for subsequent tasks.
    """
    yaml_header='fetch'
    """
    YAML header for the FetchTask, used to identify the task in configuration files as part of a ``tasks`` list.
    This header is used to declare FetchTask objects in YAML task lists.
    """
    
    def do(self):
        """ Execute the fetch task. """
        self.log_message('initiated')
        source=self.specs.get('source','pdb')
        sourceID=self.specs.get('sourceID','')
        source_format=self.specs.get('source_format','pdb')
        if source_format not in ['pdb', 'cif']:
            raise ValueError(f"Unsupported source format: {source_format}. Expected 'pdb' or 'cif'.")
        if source=='pdb':
            if source_format=='pdb':
                # parse the PDB file and register it in the pipeline context
                pdb_file=PDBParser(PDBcode=sourceID).fetch()
                if pdb_file:
                    self.ctx.register(
                        key='base_coordinates',
                        value=PDBFile(path=pdb_file),
                        value_type=PDBFile,
                        produced_by=self,
                        type='initial',
                        propagate=True
                    )
                else:
                    raise ValueError(f"Could not fetch PDB file for sourceID: {sourceID}")
            else:
                cif_file=PDBParser(PDBcode=sourceID, input_format='mmCIF').fetch()
                if cif_file:
                    self.ctx.register(
                        key='base_coordinates',
                        value=CIFFile(path=cif_file),
                        value_type=CIFFile,
                        produced_by=self,
                        type='initial',
                        propagate=True
                    )
                else:
                    raise ValueError(f"Could not fetch CIF file for sourceID: {sourceID}")
                self.log_message(f'Parsed PDB file: {sourceID}')
        elif source=='alphafold':
            pdb_file=PDBParser(alphafold=sourceID).fetch()
            if pdb_file:
                self.ctx.register(
                    key='base_coordinates',
                    value=PDBFile(path=pdb_file),
                    value_type=PDBFile,
                    produced_by=self,
                    type='initial',
                    propagate=True
                )
            else:
                raise ValueError(f"Could not fetch AlphaFold PDB file for sourceID: {sourceID}")
        elif source=='local': # expect <sourceID>.pdb or <sourceID>.cif
            local_pdb=f'{sourceID}.pdb'
            local_cif=f'{sourceID}.cif'
            if source_format=='pdb' and os.path.exists(local_pdb):
                self.ctx.register(
                    key='base_coordinates',
                    value=PDBFile(path=local_pdb),
                    value_type=PDBFile,
                    produced_by=self,
                    type='initial',
                    propagate=True
                )
            elif os.path.exists(local_cif):
                self.ctx.register(
                    key='base_cif',
                    value=CIFFile(path=local_cif),
                    value_type=CIFFile,
                    produced_by=self,
                    type='initial',
                    propagate=True
                )
            else:
                raise FileNotFoundError(f"Neither {local_pdb} nor {local_cif} found.")
        else:
            raise ValueError(f"Unsupported source type: {source}")
        self.log_message('complete')
        return self.result