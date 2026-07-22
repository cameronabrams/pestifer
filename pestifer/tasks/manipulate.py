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

    # coord-mod objtypes done in pure numpy (Phase 2); the rest (orient, crotations, irotations)
    # remain on the VMD path.
    _PYTHON_OBJTYPES = ('transrot', 'align', 'transfer_coords')

    def coormods(self):
        """
        Perform coordinate modifications based on the specifications provided in the task.
        This method retrieves the coordinate modifications from the object manager and applies each
        objtype in turn, each as an isolated load / transform / write-pdb pass that updates the
        pipeline state.  ``transrot``/``align``/``transfer_coords`` are applied in pure numpy
        (:mod:`pestifer.molecule.coordmanip`); ``orient``/``crotations``/``irotations`` still emit
        and run a VMD script.  Returns 0 on success or a non-zero error code on failure.
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
            state: StateArtifacts = self.get_current_artifact('state')
            if objtype in self._PYTHON_OBJTYPES:
                result = self._apply_python_coormods(objtype, objlist, state)
                if result != 0:
                    return result
                self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf),
                              key='state', artifact_type=StateArtifacts)
            else:
                result = self._apply_vmd_coormods(objtype, objlist, state)
                if result != 0:
                    return result
                self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf),
                              key='state', artifact_type=StateArtifacts)
                self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
                self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)
        return 0

    def _apply_python_coormods(self, objtype, objlist, state: StateArtifacts) -> int:
        """Apply a numpy-native coord-mod objtype to the current state, writing ``{basename}.pdb``."""
        from ..molecule.coordmanip import CoordManipulator, CoordManipulateError
        cm = CoordManipulator(state.psf.name, state.pdb.name)
        try:
            for obj in objlist:
                if objtype == 'transrot':
                    cm.apply_rottrans(obj)
                elif objtype == 'align':
                    cm.apply_align(obj, self._load_align_reference(obj))
                elif objtype == 'transfer_coords':
                    donor = CoordManipulator(obj.donor_psf, obj.donor_pdb)
                    cm.apply_transfer_coords(obj, donor)
        except CoordManipulateError as e:
            logger.error(f'manipulate {objtype}: {e}')
            return 1
        cm.write_pdb(f'{self.basename}.pdb')
        return 0

    def _load_align_reference(self, align: Align):
        """Resolve an Align reference (downloading from RCSB when given as a ``ref_sourceID``) and
        load it into a :class:`~pestifer.molecule.coordmanip.CoordManipulator`."""
        from ..molecule.coordmanip import CoordManipulator
        if align.ref_sourceID is not None:
            if not Path(align.effective_ref_pdb).exists():
                logger.debug(f'align: downloading reference {align.ref_sourceID} from RCSB')
                PDBParser(source_db='rcsb', source_id=align.ref_sourceID).fetch()
            self.register(align.ref_sourceID, key=f'ref_pdb_{align.ref_sourceID}', artifact_type=PDBFileArtifact)
        return CoordManipulator(align.ref_psf, align.effective_ref_pdb)

    def _apply_vmd_coormods(self, objtype, objlist, state: StateArtifacts) -> int:
        """Apply a VMD-native coord-mod objtype (orient / crotations / irotations)."""
        vm: VMDScripter = self.get_scripter('vmd')
        vm.newscript(self.basename, packages=['Orient'])
        vm.load_psf_pdb(state.psf.name, state.pdb.name, new_molid_varname='mCM')
        for obj in objlist:
            if objtype in ('irotations', 'crotations'):
                vm.write_crot(obj, molid='mCM')
            elif objtype == 'orient':
                vm.write_orient(obj, molid='mCM')
        vm.write_pdb(self.basename, 'mCM')
        vm.writescript()
        result = vm.runscript()
        # VMD's Tcl `exit N` always returns process code 0, so a guard that calls `exit 1` cannot
        # fail the task through the return code.  Detect hard-errors by scanning the log for the
        # PESTIFER-ERROR marker.
        if result == 0:
            log_path = Path(f'{self.basename}.log')
            if log_path.exists() and 'PESTIFER-ERROR' in log_path.read_text():
                logger.error(f'VMD script {self.basename}.tcl reported a PESTIFER-ERROR '
                             f'(see {log_path}); aborting manipulate task')
                result = 1
        if result != 0:
            logger.error(f'Error running VMD script {self.basename}.tcl: result {result}')
        return result

