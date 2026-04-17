# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Integration tests for DesolvateTask  (require VMD + catdcd; run with --runslow).

Fixtures are extracted on-demand from the BPTI build tarball in scratch/builds/1/.
The solvated PSF carries three DISU patch remarks; the test verifies that the
desolvated (dry) PSF preserves them.
"""
import os
import shutil
import tarfile
import unittest
from pathlib import Path

import pytest

from pestifer.core.config import Config
from pestifer.core.controller import Controller


_TARBALL = (Path(__file__).parents[3] / 'scratch' / 'builds' / '1' / 'artifacts.tar.gz').resolve()
_ARCHIVE_PREFIX = 'my_archive'
_SOLVATED_PSF = '00-04-00_solvate.psf'
_SOLVATED_PDB = '00-04-00_solvate.pdb'
_MINIMIZE_DCD = '00-05-00_md-minimize.dcd'
_DISU_REMARKS = {
    'patch DISU A:5 A:55',
    'patch DISU A:14 A:38',
    'patch DISU A:30 A:51',
}


def _extract_from_tarball(*member_basenames: str):
    """Extract named files from the build tarball into the current directory."""
    with tarfile.open(_TARBALL) as tar:
        for basename in member_basenames:
            member = tar.getmember(f'{_ARCHIVE_PREFIX}/{basename}')
            member.name = basename
            tar.extract(member, path='.')


def _psf_patch_remarks(psf_path: str) -> set[str]:
    """Return normalised patch remark texts from a PSF header."""
    remarks = set()
    with open(psf_path) as fh:
        for line in fh:
            line = line.strip()
            if '!NATOM' in line:
                break
            if line.startswith('REMARKS patch'):
                remarks.add(' '.join(line[len('REMARKS '):].split()))
    return remarks


class TestDesolvateTask(unittest.TestCase):

    def setUp(self):
        self.controller = Controller().configure(Config().configure_new(), terminate=False)

    def _setup_subdir(self, name: str):
        if os.path.exists(name):
            shutil.rmtree(name)
        os.makedirs(name)
        os.chdir(name)
        _extract_from_tarball(_SOLVATED_PSF, _SOLVATED_PDB, _MINIMIZE_DCD)

    @pytest.mark.slow
    def test_desolvate_preserves_patch_remarks(self):
        """Desolvating a solvated BPTI system must preserve DISU patch remarks."""
        self._setup_subdir('__test_desolvate_patch_remarks')
        task_list = [{'desolvate': {
            'psf': _SOLVATED_PSF,
            'pdb': _SOLVATED_PDB,
            'keepatselstr': 'protein',
            'idx_outfile': 'dry.idx',
            'psf_outfile': 'dry.psf',
            'dcd_outfile': 'dry.dcd',
            'dcd_infiles': [_MINIMIZE_DCD],
            'dcd_stride': 1,
        }}]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        self.assertTrue(Path('dry.psf').exists(), 'dry.psf was not produced')
        dry_remarks = _psf_patch_remarks('dry.psf')
        for r in _DISU_REMARKS:
            self.assertIn(r, dry_remarks, f'patch remark {r!r} missing from dry PSF')
