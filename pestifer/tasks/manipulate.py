# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`ManipulateTask` class for performing coordinate manipulations.
This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and is used for manipulating molecular coordinates.

Usage is described in the :ref:`config_ref tasks manipulate` documentation.
"""
import logging

from ..core.basetask import BaseTask
from ..core.objmanager import ObjManager

logger=logging.getLogger(__name__)

class ManipulateTask(BaseTask):
    """
    ManipulateTask class for performing coordinate manipulations.
    """
    yaml_header='manipulate'
    """
    YAML header for the ManipulateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """
    def do(self):
        """
        Execute the manipulate task.
        """
        self.log_message('initiated')
        if self.prior:
            logger.debug(f'Task {self.taskname} prior {self.prior.taskname}')
            self.inherit_state()
        logger.debug(f'manipulate {self.specs["mods"]}')
        self.objmanager=ObjManager()
        self.objmanager.ingest(self.specs['mods'])
        self.result=self.coormods()
        self.log_message('complete')
        return super().do()

    def coormods(self):
        """
        Perform coordinate modifications based on the specifications provided in the task.
        This method retrieves the coordinate modifications from the object manager,
        iterates through the specified object types and their corresponding lists,
        and applies the modifications using the VMD scripter.
        It writes the modified coordinates to a PDB file and saves the state with the updated coordinates.
        The method returns 0 on success or a non-zero error code on failure.
        """
        coord=self.objmanager.get('coord',{})
        logger.debug(f'coord {coord}')
        if not coord:
            logger.debug(f'no coormods to perform')
            return 0
        logger.debug(f'performing coormods')
        for objtype,objlist in coord.items():
            self.next_basename(objtype)
            vm=self.scripters['vmd']
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
