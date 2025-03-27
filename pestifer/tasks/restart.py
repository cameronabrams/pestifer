# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import shutil

from ..basetask import BaseTask

logger=logging.getLogger(__name__)
        
class RestartTask(BaseTask):
    """ This task only resets the workflow chain to named values of the psf, pdb, xsc, and coor files """
    yaml_header='restart'
    def do(self):
        self.log_message('initiated')
        self.next_basename()
        exts_actual=[]
        self.keepfiles=[]
        for ext in ['psf','coor','pdb','xsc']:
            fname=self.specs.get(ext,'')
            if fname:
                shutil.copy(fname,f'{self.basename}.{ext}')
                exts_actual.append(ext)
                self.keepfiles.append(fname)
        self.save_state(exts=exts_actual)
        self.log_message('complete')
        self.result=0
        return super().do()