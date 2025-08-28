# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`ManipulateTask` class for performing coordinate manipulations.
This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and is used for manipulating molecular coordinates.

Usage is described in the :ref:`config_ref tasks manipulate` documentation.
"""
import logging

from pathlib import Path

from .basetask import BaseTask
from ..scripters import VMDScripter
from ..core.objmanager import ObjManager
from ..core.artifacts import PDBFileArtifact, VMDScriptArtifact, VMDLogFileArtifact, StateArtifacts

logger = logging.getLogger(__name__)

class ManipulateTask(BaseTask):
    """
    ManipulateTask class for performing coordinate manipulations.
    """
    _yaml_header = 'manipulate'
    """
    YAML header for the ManipulateTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """
    def provision(self, packet = None):
        super().provision(packet)
        self.objmanager = ObjManager()
        self.objmanager.ingest(self.specs['mods'])

    def do(self) -> int:
        """
        Execute the manipulate task.
        """
        self.result = self.coormods()
        return self.result

    def coormods(self):
        """
        Perform coordinate modifications based on the specifications provided in the task.
        This method retrieves the coordinate modifications from the object manager,
        iterates through the specified object types and their corresponding lists,
        and applies the modifications using the VMD scripter.
        It writes the modified coordinates to a PDB file and saves the state with the updated coordinates.
        The method returns 0 on success or a non-zero error code on failure.
        """
        coord: dict = self.objmanager.get('coord', {})
        if not coord:
            logger.debug(f'no coormods to perform')
            return 0
        logger.debug(f'performing coormods')
        for objtype, objlist in coord.items():
            self.next_basename(objtype)
            vm: VMDScripter = self.get_scripter('vmd')
            vm.newscript(self.basename, packages=['Orient'])
            state: StateArtifacts = self.get_current_artifact('state')
            vm.load_psf_pdb(state.psf.name, state.pdb.name, new_molid_varname='mCM')
            for obj in objlist:
                if objtype == 'crot':
                    vm.write_crot(obj, molid='mCM')
                elif objtype == 'orient':
                    vm.write_orient(obj, molid='mCM')
                elif objtype == 'rottrans':
                    vm.write_rottrans(obj, molid='mCM')
            vm.write_pdb(self.basename,'mCM')
            vm.writescript()
            result = vm.runscript()
            if result != 0:
                return result
            self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf), key='state', artifact_type=StateArtifacts)
            # for at in [VMDScriptArtifact, VMDLogFileArtifact]:
            #     self.register(at(self.basename))
            self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
            self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)
        return 0

