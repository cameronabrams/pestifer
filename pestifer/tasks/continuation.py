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
        # self.register(PSFFileArtifact(psf.replace('.psf','')))
        # must have either coor or pdb
        coor: str = self.specs.get('coor','')
        pdb: str = self.specs.get('pdb','')
        xsc: str = self.specs.get('xsc','')
        vel: str = self.specs.get('vel','')
        if not psf:
            raise ValueError('psf file must be specified in the continuation task')
        if not coor and not pdb:
            raise ValueError('Either coor or pdb file must be specified in the continuation task')
        if coor and not pdb: 
            pdb = self.coor_to_pdb()
        if pdb and not coor: 
            coor = self.pdb_to_coor()
        self.register(NAMDMolInputDict(psf=PSFFileArtifact(psf),
                                       pdb=PDBFileArtifact(pdb),
                                       coor=NAMDCoorFileArtifact(coor),
                                       xsc=NAMDXscFileArtifact(xsc) if xsc else None,
                                       vel=NAMDVelFileArtifact(vel) if vel else None))
        self.result = 0
        return self.result