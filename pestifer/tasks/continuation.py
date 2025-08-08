# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`ContinuationTask` class for resetting the task chain to named values of the psf, pdb, xsc, and coor files.

Usage is described in the :ref:`config_ref tasks restart` documentation.
"""
import logging

from .basetask import VMDTask
from ..core.artifacts import *

logger = logging.getLogger(__name__)

class ContinuationTask(VMDTask):
    """ 
    This task only resets the task chain to named values of the psf, pdb, xsc, and coor files 
    """
    _yaml_header = 'continuation'
    
    def do(self):
        self.next_basename()
        # must have psf, either coor or pdb
        psf: str = self.specs.get('psf', '')
        if not psf:
            raise ValueError('psf file must be specified in the continuation task')
        self.register(PSFFileArtifact(psf.replace('.psf','')))
        # must have either coor or pdb
        coor: str = self.specs.get('coor','')
        pdb: str = self.specs.get('pdb','')
        if not coor and not pdb:
            raise ValueError('Either coor or pdb file must be specified in the continuation task')
        if coor:
            self.register(NAMDCoorFileArtifact(coor.replace('.coor','')))
            if not pdb: self.coor_to_pdb()
        if pdb:
            self.register(PDBFileArtifact(pdb.replace('.pdb','')))
            if not coor: self.pdb_to_coor()
        # optional xsc file
        xsc: str = self.specs.get('xsc','')
        if xsc:
            self.register(NAMDXscFileArtifact(xsc.replace('.xsc','')))
        # optional vel file
        vel: str = self.specs.get('vel','')
        if vel:
            self.register(NAMDVelFileArtifact(vel.replace('.vel','')))
        self.result = 0
        return self.result