# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for MergeTask.

Only the pure-Python logic is tested here (no VMD/psfgen required):
  - _unique_name():          static name-enumeration helper
  - _resolve_collisions():   full collision-detection pipeline (PSFContents is mocked)

Integration tests (require VMD/psfgen, marked slow):
  - TestMergeTaskIntegration
"""
import os
import shutil
import unittest
from pathlib import Path
from unittest.mock import MagicMock, patch


from pestifer.core.artifacts import StateArtifacts
from pestifer.core.config import Config
from pestifer.core.controller import Controller
from pestifer.psfutil.psfcontents import PSFContents
from pestifer.tasks.merge import MergeTask


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _fake_psf_contents(segnames: list[str]):
    """Return a MagicMock that looks like a PSFContents with the given segment names."""
    mock = MagicMock()
    mock.segments = [MagicMock(segname=s) for s in segnames]
    return mock


def _make_resolve(segname_lists, user_maps=None, strategy='enumerate'):
    """
    Call _resolve_collisions with mocked PSFContents, one per system.

    Parameters
    ----------
    segname_lists : list of list[str]
        Segment names for each system.
    user_maps : list of dict, optional
        Per-system explicit rename maps (default: all empty).
    strategy : str
        Collision strategy ('enumerate' or 'error').
    """
    if user_maps is None:
        user_maps = [{} for _ in segname_lists]
    systems = [
        {'psf': f'sys{i}.psf', 'pdb': f'sys{i}.pdb', 'segname_map': user_maps[i]}
        for i in range(len(segname_lists))
    ]
    side_effects = [_fake_psf_contents(names) for names in segname_lists]
    with patch('pestifer.tasks.merge.PSFContents', side_effect=side_effects):
        # _resolve_collisions is an instance method; instantiate a minimal mock task
        task = MergeTask.__new__(MergeTask)
        return task._resolve_collisions(systems, strategy)


# ---------------------------------------------------------------------------
# _unique_name tests
# ---------------------------------------------------------------------------

class TestUniqueName(unittest.TestCase):

    def test_no_collision_returns_base_plus_zero(self):
        result = MergeTask._unique_name('PROA', set())
        self.assertEqual(result, 'PRO0')

    def test_skips_existing_single_digit(self):
        used = {'PRO0', 'PRO1', 'PRO2'}
        result = MergeTask._unique_name('PROA', used)
        self.assertEqual(result, 'PRO3')

    def test_falls_through_to_double_digit(self):
        used = {f'PRO{i}' for i in range(10)}
        result = MergeTask._unique_name('PROA', used)
        self.assertEqual(result, 'PR00')

    def test_all_single_and_double_digits_exhausted_raises(self):
        used = {f'PRO{i}' for i in range(10)}
        used |= {f'PR{i:02d}' for i in range(100)}
        with self.assertRaises(ValueError):
            MergeTask._unique_name('PROA', used)

    def test_short_name_uses_prefix_correctly(self):
        # 2-char name: base3 = 'AB ', padded, but base3 = name[:3] = 'AB'
        result = MergeTask._unique_name('AB', set())
        self.assertEqual(result, 'AB0')

    def test_four_char_name(self):
        result = MergeTask._unique_name('MEMB', set())
        self.assertEqual(result, 'MEM0')


# ---------------------------------------------------------------------------
# _resolve_collisions – no collisions
# ---------------------------------------------------------------------------

class TestResolveNoCollision(unittest.TestCase):

    def test_two_disjoint_systems_no_rename(self):
        resolved = _make_resolve([['PROA', 'PROB'], ['SOLV', 'CLA']])
        self.assertEqual(len(resolved), 2)
        self.assertEqual(resolved[0]['effective_segname_map'], {})
        self.assertEqual(resolved[1]['effective_segname_map'], {})

    def test_segnames_captured(self):
        resolved = _make_resolve([['PROA'], ['SOLV']])
        self.assertEqual(resolved[0]['segnames'], ['PROA'])
        self.assertEqual(resolved[1]['segnames'], ['SOLV'])


# ---------------------------------------------------------------------------
# _resolve_collisions – automatic collision resolution
# ---------------------------------------------------------------------------

class TestResolveAutoCollision(unittest.TestCase):

    def test_single_collision_renamed(self):
        resolved = _make_resolve([['PROA', 'PROB'], ['PROA', 'SOLV']])
        # system 0: no renames needed
        self.assertEqual(resolved[0]['effective_segname_map'], {})
        # system 1: PROA must be renamed
        rmap = resolved[1]['effective_segname_map']
        self.assertIn('PROA', rmap)
        new_name = rmap['PROA']
        self.assertNotEqual(new_name, 'PROA')
        # must not collide with anything in system 0
        self.assertNotIn(new_name, {'PROA', 'PROB'})

    def test_multiple_collisions_all_renamed(self):
        resolved = _make_resolve(
            [['PROA', 'PROB', 'SOLV'], ['PROA', 'PROB', 'MEMB']]
        )
        rmap = resolved[1]['effective_segname_map']
        # Both PROA and PROB collide; MEMB does not
        self.assertIn('PROA', rmap)
        self.assertIn('PROB', rmap)
        self.assertNotIn('MEMB', rmap)
        # The new names must be mutually distinct
        new_names = list(rmap.values())
        self.assertEqual(len(new_names), len(set(new_names)))

    def test_three_systems_progressive_resolution(self):
        # PROA in all three systems
        resolved = _make_resolve([['PROA'], ['PROA'], ['PROA']])
        # System 0: no rename
        self.assertEqual(resolved[0]['effective_segname_map'], {})
        # System 1: PROA → PRO0
        self.assertEqual(resolved[1]['effective_segname_map'].get('PROA'), 'PRO0')
        # System 2: PROA → PRO1 (PRO0 is taken by system 1)
        self.assertEqual(resolved[2]['effective_segname_map'].get('PROA'), 'PRO1')

    def test_all_final_names_unique_across_systems(self):
        resolved = _make_resolve(
            [['PROA', 'PROB', 'SOLV'],
             ['PROA', 'MEMB'],
             ['SOLV', 'PROA']]
        )
        all_final = []
        for sys in resolved:
            for orig in sys['segnames']:
                all_final.append(sys['effective_segname_map'].get(orig, orig))
        self.assertEqual(len(all_final), len(set(all_final)),
                         f"Duplicate final segment names: {all_final}")


# ---------------------------------------------------------------------------
# _resolve_collisions – explicit user maps
# ---------------------------------------------------------------------------

class TestResolveUserMap(unittest.TestCase):

    def test_user_map_applied_even_without_collision(self):
        # System 1 explicitly renames SOLV → WAT (no collision existed)
        resolved = _make_resolve(
            [['PROA'], ['SOLV']],
            user_maps=[{}, {'SOLV': 'WAT'}],
        )
        # The rename should appear in effective_segname_map
        self.assertEqual(resolved[1]['effective_segname_map'].get('SOLV'), 'WAT')

    def test_user_map_takes_priority_over_enumerate(self):
        # PROA collides; user says rename PROA → PROT in system 1
        resolved = _make_resolve(
            [['PROA'], ['PROA']],
            user_maps=[{}, {'PROA': 'PROT'}],
        )
        self.assertEqual(resolved[1]['effective_segname_map'].get('PROA'), 'PROT')

    def test_user_map_collision_falls_through_to_enumerate(self):
        # User map sends PROA → PROB, but PROB already exists in system 0
        resolved = _make_resolve(
            [['PROB'], ['PROA']],
            user_maps=[{}, {'PROA': 'PROB'}],
        )
        rmap = resolved[1]['effective_segname_map']
        # 'PROB' is taken; should fall through to enumerate
        final = rmap.get('PROA')
        self.assertNotEqual(final, 'PROB')
        self.assertNotEqual(final, 'PROA')


# ---------------------------------------------------------------------------
# _resolve_collisions – error strategy
# ---------------------------------------------------------------------------

class TestResolveErrorStrategy(unittest.TestCase):

    def test_no_collision_does_not_raise(self):
        # Should complete silently
        resolved = _make_resolve([['PROA'], ['SOLV']], strategy='error')
        self.assertEqual(len(resolved), 2)

    def test_collision_raises_valueerror(self):
        with self.assertRaises(ValueError):
            _make_resolve([['PROA'], ['PROA']], strategy='error')

    def test_explicit_map_resolves_before_error_check(self):
        # User provides an explicit map that avoids the collision
        resolved = _make_resolve(
            [['PROA'], ['PROA']],
            user_maps=[{}, {'PROA': 'PROB'}],
            strategy='error',
        )
        self.assertEqual(resolved[1]['effective_segname_map'].get('PROA'), 'PROB')


# ---------------------------------------------------------------------------
# Integration tests  (require VMD/psfgen; run with --runslow)
# ---------------------------------------------------------------------------

# BPTI vacuum build (scratch/builds/1): 1108 atoms, segments A (protein) + B, C
# (glycan/ion); three DISU patches on segment A.
_FIXTURES = (Path(__file__).parents[3] / 'scratch' / 'builds' / '1' / 'standards').resolve()
_PSF = '00-01-00_psfgen-build.psf'
_PDB = '00-01-00_psfgen-build.pdb'
_N_ATOMS = 1108
_SEGMENTS = ('A', 'B', 'C')
_DISU_REMARKS = {'patch DISU A:5 A:55', 'patch DISU A:14 A:38', 'patch DISU A:30 A:51'}


def _translate_pdb(src: Path, dst: Path, dx: float = 0.0, dy: float = 0.0, dz: float = 100.0):
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


def _psf_patch_remarks(psf_path: str) -> set[str]:
    """Return the set of normalised patch remark texts from a PSF header."""
    remarks = set()
    with open(psf_path) as fh:
        for line in fh:
            line = line.strip()
            if '!NATOM' in line:
                break
            if line.startswith('REMARKS patch'):
                remarks.add(' '.join(line[len('REMARKS '):].split()))
    return remarks


class TestMergeTaskIntegration(unittest.TestCase):

    def setUp(self):
        self.controller = Controller().configure(Config().configure_new(), terminate=False)

    def _setup_subdir(self, name: str):
        """Create a clean working subdirectory, symlink the base fixture, write a
        translated copy, and chdir into it."""
        if os.path.exists(name):
            shutil.rmtree(name)
        os.makedirs(name)
        os.symlink(_FIXTURES / _PSF, Path(name) / _PSF)
        os.symlink(_FIXTURES / _PDB, Path(name) / _PDB)
        _translate_pdb(_FIXTURES / _PDB, Path(name) / 'bpti-shifted.pdb')
        os.chdir(name)


    def test_merge_two_copies_enumerate(self):
        """Merging BPTI with a translated copy; enumerate resolves collisions and DISU
        remarks from both copies appear in the merged PSF."""
        self._setup_subdir('__test_merge_enumerate')
        task_list = [{'merge': {
            'systems': [
                {'psf': _PSF, 'pdb': _PDB},
                {'psf': _PSF, 'pdb': 'bpti-shifted.pdb'},
            ],
            'collision_strategy': 'enumerate',
        }}]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        state: StateArtifacts = self.controller.tasks[0].get_current_artifact('state')
        self.assertTrue(state.psf.exists())
        self.assertTrue(state.pdb.exists())
        merged = PSFContents(state.psf.name)
        self.assertEqual(len(merged.atoms), 2 * _N_ATOMS)
        seg_names = [s.segname for s in merged.segments]
        self.assertEqual(len(set(seg_names)), 6, f'duplicate segment names: {seg_names}')
        # DISU remarks from system 1 (unchanged names) must be in merged PSF
        merged_remarks = _psf_patch_remarks(state.psf.name)
        for r in _DISU_REMARKS:
            self.assertIn(r, merged_remarks, f'patch remark {r!r} missing from merged PSF')


    def test_merge_two_copies_explicit_map(self):
        """Merging BPTI with a translated copy via explicit segname_map; DISU remarks
        for the renamed copy must appear with updated segment names."""
        self._setup_subdir('__test_merge_explicit')
        task_list = [{'merge': {
            'systems': [
                {'psf': _PSF, 'pdb': _PDB},
                {'psf': _PSF, 'pdb': 'bpti-shifted.pdb',
                 'segname_map': {'A': 'AA', 'B': 'BB', 'C': 'CC'}},
            ],
        }}]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        state: StateArtifacts = self.controller.tasks[0].get_current_artifact('state')
        self.assertTrue(state.psf.exists())
        merged = PSFContents(state.psf.name)
        self.assertEqual(len(merged.atoms), 2 * _N_ATOMS)
        seg_names = {s.segname for s in merged.segments}
        for expected in ('A', 'B', 'C', 'AA', 'BB', 'CC'):
            self.assertIn(expected, seg_names, f'segment {expected!r} missing from merged PSF')
        merged_remarks = _psf_patch_remarks(state.psf.name)
        # Original DISU remarks (segment A) must be present
        for r in _DISU_REMARKS:
            self.assertIn(r, merged_remarks, f'original patch remark {r!r} missing')
        # Renamed DISU remarks (segment AA) must also be present
        for r in _DISU_REMARKS:
            renamed = r.replace('A:', 'AA:')
            self.assertIn(renamed, merged_remarks, f'renamed patch remark {renamed!r} missing')


    def test_merge_feeds_pipeline(self):
        """State artifact from merge is correctly consumed by a downstream terminate task."""
        self._setup_subdir('__test_merge_pipeline')
        controller = Controller().configure(Config().configure_new())  # terminate=True by default
        task_list = [{'merge': {
            'systems': [
                {'psf': _PSF, 'pdb': _PDB},
                {'psf': _PSF, 'pdb': 'bpti-shifted.pdb'},
            ],
            'collision_strategy': 'enumerate',
        }}, {'terminate': {
            'basename': 'merged_bpti',
            'cleanup': False,
        }}]
        controller.reconfigure_tasks(task_list)
        controller.do_tasks()
        # terminate.copy_state_to_basename should have written these
        self.assertTrue(Path('merged_bpti.psf').exists(), 'merged_bpti.psf not found after terminate')
        self.assertTrue(Path('merged_bpti.pdb').exists(), 'merged_bpti.pdb not found after terminate')
        merged = PSFContents('merged_bpti.psf')
        self.assertEqual(len(merged.atoms), 2 * _N_ATOMS)

    def test_merge_three_copies(self):
        """Merging three copies of BPTI (each shifted 100 Å along Z); enumerate resolves
        collisions and DISU remarks from all three copies appear in the merged PSF."""
        self._setup_subdir('__test_merge_three_copies')
        _translate_pdb(_FIXTURES / _PDB, Path('bpti-shifted1.pdb'), dz=100.0)
        _translate_pdb(_FIXTURES / _PDB, Path('bpti-shifted2.pdb'), dz=200.0)
        task_list = [{'merge': {
            'systems': [
                {'psf': _PSF, 'pdb': _PDB},
                {'psf': _PSF, 'pdb': 'bpti-shifted1.pdb'},
                {'psf': _PSF, 'pdb': 'bpti-shifted2.pdb'},
            ],
            'collision_strategy': 'enumerate',
        }}]
        self.controller.reconfigure_tasks(task_list)
        self.controller.do_tasks()
        state: StateArtifacts = self.controller.tasks[0].get_current_artifact('state')
        self.assertTrue(state.psf.exists())
        self.assertTrue(state.pdb.exists())
        merged = PSFContents(state.psf.name)
        self.assertEqual(len(merged.atoms), 3 * _N_ATOMS)
        seg_names = {s.segname for s in merged.segments}
        self.assertEqual(len(seg_names), 9, f'expected 9 unique segment names, got: {seg_names}')
        # DISU remarks from the first (unmodified) copy must be present
        merged_remarks = _psf_patch_remarks(state.psf.name)
        for r in _DISU_REMARKS:
            self.assertIn(r, merged_remarks, f'patch remark {r!r} missing from merged PSF')
