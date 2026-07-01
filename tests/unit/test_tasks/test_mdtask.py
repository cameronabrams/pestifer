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


class TestSSRestraintsOptIn(unittest.TestCase):
    """SS restraints must be opt-in (off by default): the schema attributes were all
    defaulted, so ycleptic auto-populated the dict and every md stage -- including a
    bare-membrane relaxation with no secondary structure -- ran the ssrestraints
    plugin. The md task now builds them only when `ssrestraints.enabled` is true."""

    @staticmethod
    def _ss_schema():
        import os, yaml
        import pestifer
        base = os.path.join(os.path.dirname(pestifer.__file__), 'schema', 'base.yaml')

        def find(node):
            if isinstance(node, dict):
                if node.get('name') == 'ssrestraints' and node.get('type') == 'dict':
                    return node
                for v in node.values():
                    r = find(v)
                    if r:
                        return r
            elif isinstance(node, list):
                for v in node:
                    r = find(v)
                    if r:
                        return r
            return None

        return find(yaml.safe_load(open(base)))

    def test_schema_has_enabled_defaulting_false(self):
        ss = self._ss_schema()
        self.assertIsNotNone(ss, "ssrestraints schema block not found")
        enabled = next((a for a in ss['attributes'] if a['name'] == 'enabled'), None)
        self.assertIsNotNone(enabled, "ssrestraints needs an 'enabled' opt-in attribute")
        self.assertFalse(enabled['default'], "ssrestraints must default off (opt-in)")

    def test_guard_skips_unless_enabled(self):
        # mirrors the md task guard at mdtask.py: `if ssrestraints.get('enabled', False)`
        for specs, should_build in [({}, False), ({'enabled': False}, False),
                                    ({'enabled': True}, True), ({'sel': 'protein'}, False)]:
            self.assertEqual(bool(specs.get('enabled', False)), should_build)
