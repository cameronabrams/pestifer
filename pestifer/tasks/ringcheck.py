# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`RingCheckTask` class for checking for pierced rings in a molecular structure.
This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and is used to check for pierced rings in a molecular structure.
It identifies configurations where a ring is pierced by another segment and can optionally delete the pierced segments.
Usage is described in the :ref:`config_ref tasks ring_check` documentation.
"""
import logging

from ..core.basetask import BaseTask
from ..core.artifacts import PDBFile, PSFFile, PsfgenInputScript, PsfgenLogFile
from ..psfutil.psfring import ring_check

logger=logging.getLogger(__name__)

class RingCheckTask(BaseTask):
    """
    A class for checking for pierced rings
    """
    
    _yaml_header='ring_check'
    """
    YAML header for RingCheckTask objects.
    This header is used to declare RingCheckTask objects in YAML task lists.
    """
    
    def do(self):
        psf=self.get_current_artifact_path('psf')
        pdb=self.get_current_artifact_path('pdb')
        xsc=self.get_current_artifact_path('xsc')
        cutoff=self.specs.get('cutoff',3.5)
        segtypes=self.specs.get('segtypes',['lipid'])
        delete_these=self.specs.get('delete','piercee')
        npiercings=ring_check(psf,pdb,xsc,cutoff=cutoff,segtypes=segtypes)
        if npiercings:
            ess='s' if len(npiercings)>1 else ''
            if delete_these=="none":
                logger.debug(f'No action taken regarding {len(npiercings)} pierced-ring configuration{ess}')
                for r in npiercings:
                    logger.debug(f'  Piercing of {r["piercee"]["segname"]}-{r["piercee"]["resid"]} by {r["piercer"]["segname"]}-{r["piercer"]["resid"]}')
            else:
                self.next_basename('ring_check')
                pg=self.pipeline.get_scripter('psfgen')
                pg.newscript(self.basename)
                pg.load_project(psf.name,pdb.name)
                logger.debug(f'Deleting all {delete_these}s from {len(npiercings)} pierced-ring configuration{ess}')
                for r in npiercings:
                    logger.debug(f'   Deleting segname {r[delete_these]["segname"]} residue {r[delete_these]["resid"]}')
                    pg.addline(f'delatom {r[delete_these]["segname"]} {r[delete_these]["resid"]}')
                pg.writescript(self.basename)
                pg.runscript()
                for at in [PsfgenInputScript, PSFFile, PDBFile, PsfgenLogFile]:
                    self.register_current_artifact(at(self.basename))
        self.log_message('complete')
        self.result=0
        return super().do()