# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Tests for ManipulateTask — specifically the 'align' coormod.

Unit tests (no VMD/NAMD needed):
  - TestWriteAlign: write_align() TCL generation

Integration tests (require VMD; run with --runslow):
  - TestManipulateAlignIntegration: align a displaced system back to a reference
"""
import math
import os
import shutil
import unittest
from pathlib import Path

import pytest

from pestifer.objs.align import Align
from pestifer.scripters.vmd import VMDScripter


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

class _MockScripter(VMDScripter):
    """Minimal stub that captures addline() calls without touching the filesystem."""
    def __init__(self, molid_varname='mCM'):
        self.lines = []
        self.molid_varname = molid_varname

    def addline(self, line):
        self.lines.append(line)


def _translate_pdb(src: Path, dst: Path, dx: float, dy: float, dz: float):
    """Write a copy of a PDB with all ATOM/HETATM coordinates shifted by (dx, dy, dz)."""
    with open(src) as fh:
        lines = fh.readlines()
    with open(dst, 'w') as fh:
        for line in lines:
            if line[:6] in ('ATOM  ', 'HETATM'):
                x = float(line[30:38]) + dx
                y = float(line[38:46]) + dy
                z = float(line[46:54]) + dz
                line = f'{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}'
            fh.write(line)


def _pdb_coords(path: Path) -> list[tuple[float, float, float]]:
    """Return (x, y, z) for every ATOM/HETATM record in order."""
    coords = []
    with open(path) as fh:
        for line in fh:
            if line[:6] in ('ATOM  ', 'HETATM'):
                coords.append((
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ))
    return coords


def _rmsd(coords_a, coords_b) -> float:
    assert len(coords_a) == len(coords_b), 'atom count mismatch'
    n = len(coords_a)
    total = sum(
        (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2
        for (ax, ay, az), (bx, by, bz) in zip(coords_a, coords_b)
    )
    return math.sqrt(total / n)


# ---------------------------------------------------------------------------
# Unit tests: write_align() TCL generation
# ---------------------------------------------------------------------------

class TestWriteAlign(unittest.TestCase):

    def _lines(self, **kwargs) -> list[str]:
        vm = _MockScripter()
        a = Align(**kwargs)
        vm.write_align(a)
        return vm.lines

    def test_pdb_only_load(self):
        lines = self._lines(ref_pdb='ref.pdb')
        self.assertTrue(any('mol load pdb {ref.pdb}' in l for l in lines))
        self.assertFalse(any('mol load psf' in l for l in lines))

    def test_psf_pdb_load(self):
        lines = self._lines(ref_pdb='ref.pdb', ref_psf='ref.psf')
        self.assertTrue(any('mol load psf {ref.psf} pdb {ref.pdb}' in l for l in lines))

    def test_ref_sel_defaults_to_mobile_sel(self):
        lines = self._lines(ref_pdb='ref.pdb', mobile_sel='backbone')
        fit_line = next(l for l in lines if '_ref_fit' in l and 'atomselect' in l)
        self.assertIn('backbone', fit_line)

    def test_custom_ref_sel(self):
        lines = self._lines(ref_pdb='ref.pdb', mobile_sel='backbone', ref_sel='name CA')
        mob_line = next(l for l in lines if '_mob_fit' in l and 'atomselect' in l)
        ref_line = next(l for l in lines if '_ref_fit' in l and 'atomselect' in l)
        self.assertIn('backbone', mob_line)
        self.assertIn('name CA', ref_line)

    def test_apply_to_default_all(self):
        lines = self._lines(ref_pdb='ref.pdb')
        mover_line = next(l for l in lines if '_mover' in l and 'atomselect' in l)
        self.assertIn('"all"', mover_line)

    def test_custom_apply_to(self):
        lines = self._lines(ref_pdb='ref.pdb', apply_to='protein')
        mover_line = next(l for l in lines if '_mover' in l and 'atomselect' in l)
        self.assertIn('protein', mover_line)

    def test_congruency_check_present(self):
        lines = self._lines(ref_pdb='ref.pdb')
        self.assertTrue(any('_n_mob' in l for l in lines))
        self.assertTrue(any('_n_ref' in l for l in lines))
        self.assertTrue(any('exit 1' in l for l in lines))

    def test_congruency_error_message_includes_selections(self):
        lines = self._lines(ref_pdb='ref.pdb', mobile_sel='backbone', ref_sel='name CA')
        error_line = next((l for l in lines if 'ERROR align' in l), None)
        self.assertIsNotNone(error_line)
        self.assertIn('backbone', error_line)
        self.assertIn('name CA', error_line)

    def test_ref_mol_deleted_after_fit(self):
        lines = self._lines(ref_pdb='ref.pdb')
        fit_idx = next(i for i, l in enumerate(lines) if 'measure fit' in l)
        delete_idx = next(i for i, l in enumerate(lines) if 'mol delete' in l)
        self.assertGreater(delete_idx, fit_idx)

    def test_move_applied_after_fit(self):
        lines = self._lines(ref_pdb='ref.pdb')
        fit_idx = next(i for i, l in enumerate(lines) if 'measure fit' in l)
        move_idx = next(i for i, l in enumerate(lines) if '$_mover move' in l)
        self.assertGreater(move_idx, fit_idx)

    def test_molid_varname_used(self):
        vm = _MockScripter(molid_varname='myMol')
        vm.write_align(Align(ref_pdb='ref.pdb'))
        self.assertTrue(any('$myMol' in l for l in vm.lines))


# ---------------------------------------------------------------------------
# Integration test: full VMD run  (--runslow)
# ---------------------------------------------------------------------------

_FIXTURES = (Path(__file__).parents[2] / 'test_tasks' / 'fixtures' / 'continuation_inputs').resolve()
_PSF = 'my_6pti.psf'
_PDB = 'my_6pti.pdb'
_COOR = 'my_6pti.coor'
_RMSD_THRESHOLD = 0.1   # Å; pure-translation alignment should be essentially exact


@pytest.mark.slow
class TestManipulateAlignIntegration(unittest.TestCase):

    def setUp(self):
        from pestifer.core.config import Config
        from pestifer.core.controller import Controller
        self.controller = Controller().configure(Config().configure_new(), terminate=False)

    def _setup_subdir(self, name: str, dx: float, dy: float, dz: float):
        """Create a clean working subdirectory, symlink PSF/COOR, write displaced PDB."""
        if os.path.exists(name):
            shutil.rmtree(name)
        os.makedirs(name)
        os.symlink(_FIXTURES / _PSF,  Path(name) / _PSF)
        os.symlink(_FIXTURES / _COOR, Path(name) / _COOR)
        # Reference: original coordinates
        os.symlink(_FIXTURES / _PDB,  Path(name) / 'ref.pdb')
        # Mobile: displaced by (dx, dy, dz)
        _translate_pdb(_FIXTURES / _PDB, Path(name) / _PDB, dx, dy, dz)
        os.chdir(name)

    def _run_align(self):
        task_list = [
            {'continuation': {
                'psf': _PSF,
                'pdb': _PDB,
                'coor': _COOR,
            }},
            {'manipulate': {
                'mods': {
                    'align': [{'ref_pdb': 'ref.pdb'}],
                },
            }},
        ]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        return self.controller.tasks[1]

    def test_align_pure_translation(self):
        """Aligning a system displaced by a pure translation recovers the original coordinates."""
        self._setup_subdir('__test_align_trans', dx=25.0, dy=-10.0, dz=40.0)
        task = self._run_align()
        from pestifer.core.artifacts import StateArtifacts
        state: StateArtifacts = task.get_current_artifact('state')
        self.assertTrue(state.pdb.exists(), 'output PDB not found after manipulate/align')
        aligned_coords = _pdb_coords(Path(state.pdb.name))
        ref_coords = _pdb_coords(Path('ref.pdb'))
        rmsd = _rmsd(aligned_coords, ref_coords)
        self.assertLess(rmsd, _RMSD_THRESHOLD,
                        f'RMSD {rmsd:.4f} Å after alignment exceeds threshold {_RMSD_THRESHOLD} Å')
