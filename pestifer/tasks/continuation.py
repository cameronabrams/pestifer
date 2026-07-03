# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`ContinuationTask` class for resetting the task chain to named values of the psf, pdb, xsc, and coor files.

Usage is described in the :ref:`config_ref tasks continuation` documentation.
"""
import logging
import shutil

from .psfgen import PsfgenTask
from ..core.artifacts import *
from ..core.errors import PestiferBuildError
logger = logging.getLogger(__name__)


def _ensure_in_cwd(filepath: str) -> str:
    """Return filepath unchanged if it is already in the CWD; otherwise copy it here and return the basename."""
    if not filepath:
        return filepath
    src = Path(filepath).resolve()
    if src.parent == Path.cwd():
        return filepath
    dest = Path.cwd() / src.name
    if not dest.exists() or not dest.samefile(src):
        shutil.copy2(src, dest)
        logger.info(f'Copied {src} → {dest}')
    return src.name

class ContinuationTask(PsfgenTask):
    """ 
    This task only resets the task chain to named values of the psf, pdb, xsc, and coor files 
    """
    _yaml_header = 'continuation'

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, STATE
        # loads a pre-built system from files; an origin (would supersede any prior state).
        # Note: it provides files but NOT the in-memory `base_molecule` object.
        return TaskContract(requires=(), provides=(STATE,), discards_state=True)

    def do(self):
        self.next_basename()
        # must have psf, either coor or pdb
        psf: str = self.specs.get('psf', '')
        # must have either coor or pdb
        coor: str = self.specs.get('coor','')
        pdb: str = self.specs.get('pdb','')
        xsc: str = self.specs.get('xsc','')
        vel: str = self.specs.get('vel','')
        if not psf:
            raise PestiferBuildError('psf file must be specified in the continuation task')
        if not coor and not pdb:
            raise PestiferBuildError('Either coor or pdb file must be specified in the continuation task')
        psf = _ensure_in_cwd(psf)
        coor = _ensure_in_cwd(coor)
        pdb = _ensure_in_cwd(pdb)
        xsc = _ensure_in_cwd(xsc)
        vel = _ensure_in_cwd(vel)
        if coor and not pdb:
            pdb = self.coor_to_pdb(coor, psf)
        if pdb and not coor:
            coor = self.pdb_to_coor(pdb)
        topolines = []
        with open(psf, 'r') as f:
            for line in f:
                line = line.strip()
                if line and line.startswith("REMARKS topology"):
                    topolines.append(line)
                if '!NATOM' in line:
                    break
        if not topolines:
            raise PestiferBuildError('No topology information found in PSF file')
        streamfiles = []
        for topoline in topolines:
            tokens = topoline.split()
            toponame = os.path.basename(tokens[-1])
            if toponame.endswith('.str'):
                streamfiles.append(toponame)
        # any md task will at a minimum ensure all of charmmff's top-level parameter
        # files are copied to the CWD.  Since this PSF also lists stream files
        # in its topology remarks, the continuation task will also preemptively
        # copy them to cwd and register them as artifacts.  A PSF built outside pestifer
        # can record a stream file pestifer does not ship (e.g. VMD's NAMD-specific
        # 'toppar_water_ions_namd.str'); copy_charmmfile_local only warns and writes
        # nothing in that case, so registering it anyway would later emit a
        # 'topology <missing>' line that aborts psfgen.  Keep only the ones that resolve.
        CC = self.resource_manager.charmmff_content
        available_streamfiles = []
        for streamfile in streamfiles:
            CC.copy_charmmfile_local(streamfile)
            if os.path.exists(streamfile):
                available_streamfiles.append(streamfile)
            else:
                logger.warning(
                    f'continuation: topology stream file {streamfile!r} recorded in '
                    f'{psf} is not available in pestifer\'s force field and was not found '
                    f'locally; dropping it so it cannot abort a downstream psfgen run. If '
                    f'this PSF genuinely needs it, place the file in the run directory.')

        self.register([CharmmffStreamFileArtifact(x) for x in available_streamfiles], key='charmmff_streamfiles', artifact_type=CharmmffStreamFileArtifacts)

        self.register(dict(
                        psf=PSFFileArtifact(psf),
                        pdb=PDBFileArtifact(pdb),
                        coor=NAMDCoorFileArtifact(coor),
                        xsc=NAMDXscFileArtifact(xsc) if xsc else None,
                        vel=NAMDVelFileArtifact(vel) if vel else None), 
                        key='state', artifact_type=StateArtifacts)
        self.update_molecule()
        self.result = 0
        return self.result