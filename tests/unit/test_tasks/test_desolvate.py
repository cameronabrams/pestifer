# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Integration tests for DesolvateTask  (require VMD + catdcd; run with --runslow).

Fixtures are copied from the tracked BPTI build in fixtures/bpti_build/.
The solvated PSF carries three DISU patch remarks; the test verifies that the
desolvated (dry) PSF preserves them.
"""
import os
import shutil
import unittest
from pathlib import Path

import pytest

from pestifer.core.config import Config
from pestifer.core.controller import Controller

pytestmark = pytest.mark.needs_tools


# TRACKED fixtures (see the directory's README).  These used to be extracted from a 37 MB
# `artifacts.tar.gz` under the gitignored `scratch/builds/1/`, so the module passed only on a
# machine that happened to have run that build and raised FileNotFoundError on every fresh
# clone.  CI never saw it: the module is `needs_tools` and skips there.
_FIXTURES = (Path(__file__).parent / 'fixtures' / 'bpti_build').resolve()
_SOLVATED_PSF = '00-04-00_solvate.psf'
_SOLVATED_PDB = '00-04-00_solvate.pdb'
_MINIMIZE_DCD = '00-05-00_md-minimize.dcd'
_DISU_REMARKS = {
    'patch DISU A:5 A:55',
    'patch DISU A:14 A:38',
    'patch DISU A:30 A:51',
}


def _stage_fixtures(*basenames: str):
    """Copy named fixture files into the current directory."""
    for basename in basenames:
        shutil.copy(_FIXTURES / basename, basename)


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
        _stage_fixtures(_SOLVATED_PSF, _SOLVATED_PDB, _MINIMIZE_DCD)


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
