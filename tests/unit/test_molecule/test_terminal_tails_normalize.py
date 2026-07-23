import unittest

import pestifer.molecule.asymmetricunit as au
from pestifer.molecule.asymmetricunit import normalize_terminal_tails


class TestNormalizeTerminalTails(unittest.TestCase):
    """The grouped terminal_tails block folds into the legacy flat keys (and unions with them)."""

    def setUp(self):
        # the deprecation warning is one-time via a module global; reset it per test
        au._terminal_tails_deprecation_warned = False

    def test_grouped_folds_into_legacy_keys(self):
        seq = {'terminal_tails': {'n': ['T'], 'c': ['B', 'D'], 'all': False}}
        normalize_terminal_tails(seq)
        self.assertEqual(seq['build_zero_occupancy_N_termini'], ['T'])
        self.assertEqual(seq['build_zero_occupancy_C_termini'], ['B', 'D'])
        self.assertFalse(seq['include_terminal_loops'])

    def test_all_maps_to_include_terminal_loops(self):
        seq = {'terminal_tails': {'all': True}}
        normalize_terminal_tails(seq)
        self.assertTrue(seq['include_terminal_loops'])

    def test_grouped_and_legacy_union(self):
        seq = {'terminal_tails': {'n': ['T']},
               'build_zero_occupancy_N_termini': ['X'],
               'build_zero_occupancy_C_termini': ['B']}
        normalize_terminal_tails(seq)
        self.assertEqual(seq['build_zero_occupancy_N_termini'], ['T', 'X'])   # sorted union
        self.assertEqual(seq['build_zero_occupancy_C_termini'], ['B'])

    def test_empty_yields_present_default_keys(self):
        # robustness: the flat list keys must always exist afterward (daughter-chain propagation
        # relies on it); an empty/absent declaration yields empty lists, not KeyErrors.
        seq = {}
        normalize_terminal_tails(seq)
        self.assertEqual(seq['build_zero_occupancy_N_termini'], [])
        self.assertEqual(seq['build_zero_occupancy_C_termini'], [])
        self.assertFalse(seq['include_terminal_loops'])

    def test_legacy_only_warns_once(self):
        seq = {'build_zero_occupancy_C_termini': ['B', 'D', 'F']}
        with self.assertLogs('pestifer.molecule.asymmetricunit', level='WARNING') as cm:
            normalize_terminal_tails(seq)
        self.assertIn('deprecated', ''.join(cm.output))
        self.assertEqual(seq['build_zero_occupancy_C_termini'], ['B', 'D', 'F'])

    def test_grouped_only_does_not_warn(self):
        seq = {'terminal_tails': {'c': ['B']}}
        with self.assertNoLogs('pestifer.molecule.asymmetricunit', level='WARNING'):
            normalize_terminal_tails(seq)


if __name__ == '__main__':
    unittest.main()
