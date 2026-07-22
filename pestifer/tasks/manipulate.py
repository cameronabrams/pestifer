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
from ..core.objmanager import ObjManager
from ..core.artifacts import PDBFileArtifact, StateArtifacts
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
        Every coord-mod objtype -- ``transrot``/``align``/``transfer_coords``/``orient``/
        ``irotations``/``crotations`` -- is applied in pure numpy
        (:mod:`pestifer.molecule.coordmanip`), each as an isolated load / transform / write-pdb
        pass that updates the pipeline state.  Returns 0 on success or a non-zero error code.
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
            result = self._apply_python_coormods(objtype, objlist, state)
            if result != 0:
                return result
            self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf),
                          key='state', artifact_type=StateArtifacts)
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
                    cm.apply_transfer_coords(obj, CoordManipulator(obj.donor_psf, obj.donor_pdb))
                elif objtype in ('irotations', 'crotations'):
                    cm.apply_crot(obj)
                elif objtype == 'orient':
                    cm.apply_orient(obj)
                else:
                    logger.error(f'manipulate: unknown coord-mod objtype {objtype!r}')
                    return 1
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

