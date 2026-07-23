# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Tests for the shared glycan-pendant ring-piercing resolver (pestifer.psfutil.ring_resolve),
used both by RingCheckTask and by the psfgen task at glycan graft time.

The rotation path itself (clearing a real glycan-bond piercing) is validated end-to-end by a real
glycosylated build; here we test the detection/partition wiring: only a glycan-*bond* piercer of a
protein/glycan ring is resolvable by rotating the glycan, everything else is returned unresolved.
"""
import os
import unittest
from pathlib import Path
from unittest import mock

from pestifer.psfutil import ring_resolve

# reuse the ring_check fixture 5: a glycan pyranose ring pierced by a *lipid* bond
_FIX5 = (Path(__file__).parent.parent / 'test_tasks' / 'test_ringcheck' / '5').resolve()


class TestIsGlycanPierce(unittest.TestCase):
    def _p(self, piercee_type, piercer_type, bond_serials=(1, 2)):
        return {'piercee': {'segtype': piercee_type, 'segname': 'A', 'resid': 1},
                'piercer': {'segtype': piercer_type, 'segname': 'B', 'resid': 2,
                            'bond_serials': list(bond_serials)}}

    def test_glycan_bond_through_protein_ring_is_resolvable(self):
        self.assertTrue(ring_resolve._is_glycan_pierce(self._p('protein', 'glycan')))

    def test_glycan_bond_through_glycan_ring_is_resolvable(self):
        self.assertTrue(ring_resolve._is_glycan_pierce(self._p('glycan', 'glycan')))

    def test_lipid_piercer_not_resolvable(self):
        self.assertFalse(ring_resolve._is_glycan_pierce(self._p('glycan', 'lipid')))

    def test_glycan_piercer_without_bond_serials_not_resolvable(self):
        p = self._p('protein', 'glycan')
        p['piercer']['bond_serials'] = []
        self.assertFalse(ring_resolve._is_glycan_pierce(p))


class TestResolveGlycanPiercings(unittest.TestCase):

    def test_no_piercings_returns_none(self):
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[]):
            fixed, unresolved = ring_resolve.resolve_glycan_piercings('x.psf', 'x.pdb', None, 'out')
        self.assertIsNone(fixed)
        self.assertEqual(unresolved, [])

    def test_non_glycan_piercer_left_unresolved_no_rotation(self):
        # fixture 5 is a glycan ring pierced by a LIPID bond -> the glycan-pendant path cannot
        # (and must not) touch it; it comes back unresolved with no rotation attempted
        psf, pdb, xsc = (str(_FIX5 / 'test.psf'), str(_FIX5 / 'test.pdb'), str(_FIX5 / 'test.xsc'))
        with mock.patch.object(ring_resolve, 'try_glycan_pendant') as rot:
            fixed, unresolved = ring_resolve.resolve_glycan_piercings(
                psf, pdb, xsc, 'out', segtypes=['lipid', 'glycan'])
            rot.assert_not_called()
        self.assertIsNone(fixed)
        self.assertEqual(len(unresolved), 1)
        self.assertEqual(unresolved[0]['piercer']['segtype'], 'lipid')

    def test_resolvable_piercing_is_rotated(self):
        piercing = {'piercee': {'segtype': 'protein', 'segname': 'A', 'resid': 5, 'resname': 'PRO'},
                    'piercer': {'segtype': 'glycan', 'segname': 'G', 'resid': 3, 'resname': 'BGNA',
                                'bond_serials': [10, 11]}}
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[piercing]), \
             mock.patch.object(ring_resolve, 'RingChecker') as RC, \
             mock.patch.object(ring_resolve, 'cell_from_xsc', return_value=(None, None)), \
             mock.patch.object(ring_resolve, 'try_glycan_pendant', return_value='fixed.pdb') as rot:
            RC.return_value.load_coords.return_value = mock.Mock()
            fixed, unresolved = ring_resolve.resolve_glycan_piercings('x.psf', 'x.pdb', None, 'out')
            rot.assert_called_once()
        self.assertEqual(fixed, 'fixed.pdb')
        self.assertEqual(unresolved, [])

    def test_resolvable_but_rotation_fails_is_unresolved(self):
        piercing = {'piercee': {'segtype': 'glycan', 'segname': 'A', 'resid': 5, 'resname': 'BMAN'},
                    'piercer': {'segtype': 'glycan', 'segname': 'G', 'resid': 3, 'resname': 'BGNA',
                                'bond_serials': [10, 11]}}
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[piercing]), \
             mock.patch.object(ring_resolve, 'RingChecker') as RC, \
             mock.patch.object(ring_resolve, 'cell_from_xsc', return_value=(None, None)), \
             mock.patch.object(ring_resolve, 'try_glycan_pendant', return_value=None):
            RC.return_value.load_coords.return_value = mock.Mock()
            fixed, unresolved = ring_resolve.resolve_glycan_piercings('x.psf', 'x.pdb', None, 'out')
        self.assertIsNone(fixed)
        self.assertEqual(len(unresolved), 1)


class TestPsfgenGating(unittest.TestCase):
    """PsfgenTask.resolve_glycan_piercings gating: toggle off / no glycans -> no work."""

    def _task(self, *, check=True, nglycans=1):
        from pestifer.tasks.psfgen import PsfgenTask
        t = PsfgenTask.__new__(PsfgenTask)
        t.specs = {'source': {'sequence': {'glycans': {'declash': {'check_piercings': check}}}}}
        t.declash_counts = {'glycan': nglycans}
        return t

    def test_toggle_off_skips(self):
        t = self._task(check=False)
        with mock.patch('pestifer.psfutil.ring_resolve.resolve_glycan_piercings') as m:
            t.resolve_glycan_piercings()
            m.assert_not_called()

    def test_no_glycans_skips(self):
        t = self._task(check=True, nglycans=0)
        with mock.patch('pestifer.psfutil.ring_resolve.resolve_glycan_piercings') as m:
            t.resolve_glycan_piercings()
            m.assert_not_called()


class TestPiercerWithin(unittest.TestCase):
    """_piercer_within: a piercing is mover-caused iff its bond lies wholly in the mover set."""

    def test_bond_within_movers(self):
        self.assertTrue(ring_resolve._piercer_within({'piercer': {'bond_serials': [10, 11]}}, {10, 11, 12}))

    def test_bond_partly_outside_movers(self):
        self.assertFalse(ring_resolve._piercer_within({'piercer': {'bond_serials': [10, 99]}}, {10, 11}))

    def test_no_bond_serials(self):
        self.assertFalse(ring_resolve._piercer_within({'piercer': {}}, {1, 2}))


class TestResolveModelBuiltPiercings(unittest.TestCase):

    def test_no_piercings_returns_none(self):
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[]):
            fixed, unresolved = ring_resolve.resolve_model_built_piercings(
                'x.psf', 'x.pdb', None, 'out', {1, 2})
        self.assertIsNone(fixed)
        self.assertEqual(unresolved, [])

    def test_non_mover_piercing_is_ignored(self):
        # a piercing whose bond is NOT within the mover set is another pass's concern:
        # it must NOT be rotated here and must NOT be reported as unresolved.
        p = {'piercee': {'segtype': 'protein', 'segname': 'A', 'resid': 5},
             'piercer': {'segtype': 'glycan', 'segname': 'G', 'resid': 3, 'bond_serials': [100, 101]}}
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[p]), \
             mock.patch.object(ring_resolve, 'try_glycan_pendant') as rot:
            fixed, unresolved = ring_resolve.resolve_model_built_piercings(
                'x.psf', 'x.pdb', None, 'out', {10, 11})
            rot.assert_not_called()
        self.assertIsNone(fixed)
        self.assertEqual(unresolved, [])

    def test_mover_piercing_is_rotated(self):
        # a tail backbone bond (serials within the mover set) spearing a proline ring is rotated
        p = {'piercee': {'segtype': 'protein', 'segname': 'A', 'resid': 5, 'resname': 'PRO'},
             'piercer': {'segtype': 'protein', 'segname': 'A', 'resid': 60, 'resname': 'GLY',
                         'bond_serials': [10, 11]}}
        with mock.patch.object(ring_resolve, 'ring_check', return_value=[p]), \
             mock.patch.object(ring_resolve, 'RingChecker') as RC, \
             mock.patch.object(ring_resolve, 'cell_from_xsc', return_value=(None, None)), \
             mock.patch.object(ring_resolve, 'try_glycan_pendant', return_value='fixed.pdb') as rot:
            RC.return_value.load_coords.return_value = mock.Mock()
            fixed, unresolved = ring_resolve.resolve_model_built_piercings(
                'x.psf', 'x.pdb', None, 'out', {10, 11}, mover_kind='tail')
            rot.assert_called_once()
            # the mover set is forwarded so only tail atoms may move
            self.assertEqual(rot.call_args.kwargs.get('mover_serials'), {10, 11})
        self.assertEqual(fixed, 'fixed.pdb')
        self.assertEqual(unresolved, [])


if __name__ == '__main__':
    unittest.main()
