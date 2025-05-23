# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging

from ..basetask import BaseTask
from ..objmanager import ObjManager

logger=logging.getLogger(__name__)

class ManipulateTask(BaseTask):
    yaml_header='manipulate'
    def do(self):
        self.log_message('initiated')
        if self.prior:
            logger.debug(f'Task {self.taskname} prior {self.prior.taskname}')
            self.inherit_state()
        logger.debug(f'manipulate {self.specs["mods"]}')
        self.objmanager=ObjManager()
        self.objmanager.injest(self.specs['mods'])
        self.result=self.coormods()
        self.log_message('complete')
        return super().do()

    def coormods(self):
        coord=self.objmanager.get('coord',{})
        logger.debug(f'coord {coord}')
        if not coord:
            logger.debug(f'no coormods to perform')
            return 0
        logger.debug(f'performing coormods')
        for objtype,objlist in coord.items():
            self.next_basename(objtype)
            vm=self.writers['vmd']
            vm.newscript(self.basename,packages=['Orient'])
            psf=self.statevars['psf']
            pdb=self.statevars['pdb']
            vm.load_psf_pdb(psf,pdb,new_molid_varname='mCM')
            objlist.write_TcL(vm)
            vm.write_pdb(self.basename,'mCM')
            vm.writescript()
            result=vm.runscript()
            if result!=0:
                return result
            self.save_state(exts=['pdb'])
        return 0
