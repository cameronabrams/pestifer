# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`ManipulateTask` class for performing coordinate manipulations.
This class is a descendant of the :class:`BaseTask <pestifer.tasks.basetask.BaseTask>` class and is used for manipulating molecular coordinates.

Usage is described in the :ref:`config_ref tasks manipulate` documentation.
"""
import logging

from pathlib import Path
from pidibble.pdbparse import PDBParser

from .basetask import BaseTask
from ..scripters import VMDScripter
from ..core.objmanager import ObjManager
from ..core.artifacts import PDBFileArtifact, VMDScriptArtifact, VMDLogFileArtifact, StateArtifacts
from ..objs.align import Align

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
            if not objlist:
                continue
            self.next_basename(objtype)
            vm: VMDScripter = self.get_scripter('vmd')
            vm.newscript(self.basename, packages=['Orient'])
            state: StateArtifacts = self.get_current_artifact('state')
            vm.load_psf_pdb(state.psf.name, state.pdb.name, new_molid_varname='mCM')
            for obj in objlist:
                if objtype in ('irotations', 'crotations'):
                    vm.write_crot(obj, molid='mCM')
                elif objtype == 'orient':
                    vm.write_orient(obj, molid='mCM')
                elif objtype == 'transrot':
                    vm.write_rottrans(obj, molid='mCM')
                elif objtype == 'align':
                    if isinstance(obj, Align) and obj.ref_sourceID is not None:
                        ref_pdb_path = Path(f'{obj.ref_sourceID}.pdb')
                        if not ref_pdb_path.exists():
                            logger.debug(f'align: downloading reference {obj.ref_sourceID} from RCSB')
                            PDBParser(source_db='rcsb', source_id=obj.ref_sourceID).fetch()
                        self.register(obj.ref_sourceID, key=f'ref_pdb_{obj.ref_sourceID}', artifact_type=PDBFileArtifact)
                    vm.write_align(obj)
                elif objtype == 'transfer_coords':
                    vm.write_transfer_coords(obj)
            vm.write_pdb(self.basename, 'mCM')
            vm.writescript()
            result = vm.runscript()
            # VMD's Tcl `exit N` always returns process code 0, so a guard that calls `exit 1`
            # (e.g. the transrot disconnection guard) cannot fail the task through the return
            # code.  Detect such hard-errors by scanning the log for the PESTIFER-ERROR marker.
            if result == 0:
                log_path = Path(f'{self.basename}.log')
                if log_path.exists() and 'PESTIFER-ERROR' in log_path.read_text():
                    logger.error(f'VMD script {self.basename}.tcl reported a PESTIFER-ERROR '
                                 f'(see {log_path}); aborting manipulate task')
                    result = 1
            if result != 0:
                logger.error(f'Error running VMD script {self.basename}.tcl: result {result}')
                return result
            self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf), key='state', artifact_type=StateArtifacts)
            # for at in [VMDScriptArtifact, VMDLogFileArtifact]:
            #     self.register(at(self.basename))
            self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
            self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)
        return 0

