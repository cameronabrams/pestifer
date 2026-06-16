# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import unittest

from pestifer.core.errors import PestiferBuildError
from pestifer.tasks import mdtask
from pestifer.tasks.mdtask import normalize_ensemble


class TestNormalizeEnsemble(unittest.TestCase):
    """'NPgT' is the preferred name for the constant-x:y-ratio ensemble; 'NPAT' is the
    accepted-but-deprecated alias that maps to the same behavior with a one-time warning."""

    def setUp(self):
        mdtask._npat_deprecation_warned = False   # reset the once-only warning latch

    def test_npgt_variants_canonicalize(self):
        for name in ('NPgT', 'npgt', 'NPGT', 'NpGt'):
            self.assertEqual(normalize_ensemble(name), 'npgt')

    def test_plain_ensembles(self):
        self.assertEqual(normalize_ensemble('NPT'), 'npt')
        self.assertEqual(normalize_ensemble('nvt'), 'nvt')
        self.assertEqual(normalize_ensemble('Minimize'), 'minimize')

    def test_npat_aliases_to_npgt(self):
        self.assertEqual(normalize_ensemble('NPAT'), 'npgt')   # same behavior as NPgT
        self.assertEqual(normalize_ensemble('npat'), 'npgt')

    def test_npat_warns_once(self):
        with self.assertLogs('pestifer.tasks.mdtask', level='WARNING') as cm:
            normalize_ensemble('NPAT')
        self.assertTrue(any('NPAT' in m and 'NPgT' in m for m in cm.output))
        # second use must not warn again (latch set)
        with self.assertRaises(AssertionError):
            with self.assertLogs('pestifer.tasks.mdtask', level='WARNING'):
                normalize_ensemble('NPAT')

    def test_npgt_does_not_warn(self):
        with self.assertRaises(AssertionError):
            with self.assertLogs('pestifer.tasks.mdtask', level='WARNING'):
                normalize_ensemble('NPgT')

    def test_invalid_raises(self):
        with self.assertRaises(PestiferBuildError):
            normalize_ensemble('NPGAMMAT')
        with self.assertRaises(PestiferBuildError):
            normalize_ensemble('foo')
