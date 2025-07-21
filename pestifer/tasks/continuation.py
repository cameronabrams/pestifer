# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`ContinuationTask` class for resetting the task chain to named values of the psf, pdb, xsc, and coor files.

Usage is described in the :ref:`config_ref tasks restart` documentation.
"""
import logging

from ..core.basetask import VMDTask
from ..core.pipeline import PDBFile, COORFile, VELFile, XSCFile, PSFFile
logger=logging.getLogger(__name__)

class ContinuationTask(VMDTask):
    """ 
    This task only resets the task chain to named values of the psf, pdb, xsc, and coor files 
    """
    yaml_header='continuation'
    def do(self):
        self.log_message('initiated')
        self.next_basename()
        # must have psf, either coor or pdb
        psf = self.specs.get('psf','')
        if not psf:
            raise ValueError('psf file must be specified in the continuation task')
        self.register_current_artifact('psf', PSFFile(path=psf))
        # must have either coor or pdb
        coor = self.specs.get('coor','')
        pdb = self.specs.get('pdb','')
        if not coor and not pdb:
            raise ValueError('Either coor or pdb file must be specified in the continuation task')
        if coor:
            self.register_current_artifact('coor', COORFile(path=coor))
            self.coor_to_pdb()
        if pdb:
            self.register_current_artifact('pdb', PDBFile(path=pdb))
            self.pdb_to_coor()
        # optional xsc file
        xsc = self.specs.get('xsc','')
        if xsc:
            self.register_current_artifact('xsc', XSCFile(path=xsc))
        # optional vel file
        vel = self.specs.get('vel','')
        if vel:
            self.register_current_artifact('vel', VELFile(path=vel))
        self.log_message('complete')
        self.result = 0
        return super().do()