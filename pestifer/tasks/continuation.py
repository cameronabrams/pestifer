# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`ContinuationTask` class for resetting the task chain to named values of the psf, pdb, xsc, and coor files.

Usage is described in the :ref:`config_ref tasks continuation` documentation.
"""
import logging

from .psfgen import PsfgenTask
from ..core.artifacts import *
logger = logging.getLogger(__name__)

class ContinuationTask(PsfgenTask):
    """ 
    This task only resets the task chain to named values of the psf, pdb, xsc, and coor files 
    """
    _yaml_header = 'continuation'
    
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
            raise ValueError('psf file must be specified in the continuation task')
        if not coor and not pdb:
            raise ValueError('Either coor or pdb file must be specified in the continuation task')
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
            raise ValueError('No topology information found in PSF file')
        streamfiles = []
        for topoline in topolines:
            tokens = topoline.split()
            toponame = os.path.basename(tokens[-1])
            if toponame.endswith('.str'):
                streamfiles.append(toponame)
        # any md task will at a minimum ensure all of charmmff's top-level parameter
        # files are copied to the CWD.  Since this PSF also lists stream files
        # in its topology remarks, the continuation task will also preemptively
        # copy them to cwd and register them as artifacts
        CC = self.resource_manager.charmmff_content
        for streamfile in streamfiles:
            CC.copy_charmmfile_local(streamfile)
    
        self.register([CharmmffStreamFileArtifact(x) for x in streamfiles], key='charmmff_streamfiles', artifact_type=CharmmffStreamFileArtifacts)

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