# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for MergeTask.

Only the pure-Python logic is tested here (no VMD/psfgen required):
  - _unique_name():          static name-enumeration helper
  - _resolve_collisions():   full collision-detection pipeline (PSFContents is mocked)
"""
import unittest
from unittest.mock import MagicMock, patch

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
