# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`InitiateTask` class for initiating a pestifer build.
"""

from ..core.basetask import BaseTask
from ..core.pipeline import PipelineContext, PDBFile, CIFFile
from pidibble.pdbparse import PDBParser
import os

class InitiateTask(BaseTask):
    """
    InitiateTask class for initiating a pestifer build.
    This class inherits from the :class:`BaseTask` class and is used to set up the initial state of the build.
    It prepares the environment, initializes necessary configurations, and sets up the context for subsequent tasks.
    """
    yaml_header='initiate'
    """
    YAML header for the InitiateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    This header is used to declare InitiateTask objects in YAML task lists.
    """
    
    def do(self):
        """ Execute the initiate task. """
        self.log_message('initiated')
        # verify that the pipeline context is empty; i.e., no artifacts registered
        if self.ctx.artifacts:
            raise RuntimeError("Pipeline context is not empty. Artifacts already registered.")
        source=self.specs.get('source','pdb')
        sourceID=self.specs.get('sourceID','')
        if source=='pdb':
            # parse the PDB file and register it in the pipeline context
            pdb_file=PDBParser(PDBcode=sourceID).fetch()
            if pdb_file:
                self.ctx.register(
                    key='base_pdb',
                    file_obj=PDBFile(path=pdb_file),
                    produced_by=self.__class__.__name__,
                    type='initial',
                    propagate=True
                )
            else:
                cif_file=PDBParser(PDBcode=sourceID, input_format='mmCIF').fetch()
                if cif_file:
                    self.ctx.register(
                        key='base_cif',
                        file_obj=CIFFile(path=cif_file),
                        produced_by=self.__class__.__name__,
                        type='initial',
                        propagate=True
                    )
                else:
                    raise ValueError(f"Could not fetch PDB or CIF file for sourceID: {sourceID}")
                self.log_message(f'Parsed PDB file: {sourceID}')
        elif source=='alphafold':
            pdb_file=PDBParser(alphafold=sourceID).fetch()
            if pdb_file:
                self.ctx.register(
                    key='base_pdb',
                    file_obj=PDBFile(path=pdb_file),
                    produced_by=self.__class__.__name__,
                    type='initial',
                    propagate=True
                )
            else:
                raise ValueError(f"Could not fetch AlphaFold PDB file for sourceID: {sourceID}")
        elif source=='local': # expect <sourceID>.pdb or <sourceID>.cif
            local_pdb=f'{sourceID}.pdb'
            local_cif=f'{sourceID}.cif'
            if os.path.exists(local_pdb):
                self.ctx.register(
                    key='base_pdb',
                    file_obj=PDBFile(path=local_pdb),
                    produced_by=self.__class__.__name__,
                    type='initial',
                    propagate=True
                )
            elif os.path.exists(local_cif):
                self.ctx.register(
                    key='base_cif',
                    file_obj=CIFFile(path=local_cif),
                    produced_by=self.__class__.__name__,
                    type='initial',
                    propagate=True
                )
            else:
                raise FileNotFoundError(f"Neither {local_pdb} nor {local_cif} found.")
        else:
            raise ValueError(f"Unsupported source type: {source}")
        self.log_message('complete')
        return self.result