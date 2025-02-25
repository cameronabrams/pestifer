# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging

from ..basetask import BaseTask
from ..ring import ring_check

logger=logging.getLogger(__name__)

class RingCheckTask(BaseTask):
    """ A class for checking for pierced rings
    
    Attributes
    ----------
    yaml_header(str) 

    Methods
    -------
    do(): 
        Based on specs, executes ring piercing check and deletes offending residues

    """
    yaml_header='ring_check'

    def do(self):
        self.log_message('initiated')
        self.inherit_state()
        psf=self.statevars.get('psf',None)
        pdb=self.statevars.get('pdb',None)
        xsc=self.statevars.get('xsc',None)
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
                pg=self.writers['psfgen']
                pg.newscript(self.basename)
                pg.load_project(psf,pdb)
                logger.debug(f'Deleting all {delete_these}s from {len(npiercings)} pierced-ring configuration{ess}')
                for r in npiercings:
                    logger.debug(f'   Deleting segname {r[delete_these]["segname"]} residue {r[delete_these]["resid"]}')
                    pg.addline(f'delatom {r[delete_these]["segname"]} {r[delete_these]["resid"]}')
                pg.writescript(self.basename)
                pg.runscript()
                self.save_state(exts=['psf','pdb'])
        self.log_message('complete')
        self.result=0
        return super().do()