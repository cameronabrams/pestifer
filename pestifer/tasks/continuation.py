# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`ContinuationTask` class for resetting the task chain to named values of the psf, pdb, xsc, and coor files.

Usage is described in the :ref:`config_ref tasks restart` documentation.
"""
import logging

from ..core.basetask import BaseTask
from ..core.pipeline import PDBFile, COORFile, VELFile, XSCFile, PSFFile
logger=logging.getLogger(__name__)

class ContinuationTask(BaseTask):
    """ 
    This task only resets the task chain to named values of the psf, pdb, xsc, and coor files 
    """
    yaml_header='continuation'
    def do(self):
        self.log_message('initiated')
        self.next_basename()
        for ext,objtype in zip(['psf','coor','pdb','xsc','vel'],[PSFFile,COORFile,PDBFile,XSCFile,VELFile]):
            fname=self.specs.get(ext,'')
            if fname:
                self.ctx.register(
                    key=ext,
                    value=objtype(path=fname),
                    value_type=objtype,
                    produced_by=self,
                    type='initial',
                    propagate=True
                )
        self.log_message('complete')
        self.result=0
        return super().do()