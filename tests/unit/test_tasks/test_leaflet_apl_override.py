# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A supplied per-leaflet APL replaces the calibrated one, and does not replace the calibration.

``make_membrane_system`` derives each leaflet's preferred area per lipid from a relaxed symmetric
calibration patch.  A small patch systematically mis-states the APL of some mixtures: measured on
a PIP2 mix, 45.63 A^2 on the patch against 43.47 A^2 on a 700-lipid bilayer at 310 K, a ~10%
error that left one leaflet ~5% under-filled (674 placed against ~708).  Someone who has already
equilibrated a large bilayer of the same composition holds a better number than the task can
derive, so ``bilayer.leaflet_apl`` lets them supply it.

The design constraint is the reason this is a small change rather than a bypass: **pestifer
cannot tell whether a supplied number came from a converged run.**  An override therefore makes
an unconverged APL easier to use, not more correct -- and the failure it is most often reached
for is a calibration cell that had not equilibrated, which the task's existing area-drift warning
already names.  So the calibration still runs, its convergence is still checked, and the log
carries the supplied value, the value it displaced, and the percentage between them.
"""
import inspect
import unittest

from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask

_ov = MakeMembraneSystemTask._apply_leaflet_apl_override


class TestTheOverride(unittest.TestCase):

    def test_no_spec_leaves_both_calibrated(self):
        for bs in ({}, None, {'leaflet_apl': {}}, {'leaflet_apl': None}):
            self.assertEqual(_ov(bs, 41.58, 43.47), (41.58, 43.47), f'bs={bs!r}')

    def test_both_leaflets_can_be_supplied(self):
        self.assertEqual(_ov({'leaflet_apl': {'upper': 41.58, 'lower': 43.47}}, 41.32, 45.63),
                         (41.58, 43.47))

    def test_one_leaflet_can_be_supplied_and_the_other_stays_calibrated(self):
        """The asymmetric case is the point: her raft leaflet agreed, the PIP2 one did not."""
        self.assertEqual(_ov({'leaflet_apl': {'lower': 43.47}}, 41.32, 45.63), (41.32, 43.47))
        self.assertEqual(_ov({'leaflet_apl': {'upper': 41.58}}, 41.32, 45.63), (41.58, 45.63))

    def test_a_supplied_value_is_coerced_to_float(self):
        up, lo = _ov({'leaflet_apl': {'upper': '41.58', 'lower': 43}}, 41.32, 45.63)
        self.assertIsInstance(up, float)
        self.assertIsInstance(lo, float)
        self.assertAlmostEqual(up, 41.58)
        self.assertAlmostEqual(lo, 43.0)

    def test_a_nonsensical_APL_is_rejected_rather_than_built(self):
        for bad in (0, -1.0, 0.0):
            with self.assertRaises(ValueError):
                _ov({'leaflet_apl': {'upper': bad}}, 41.32, 45.63)

    def test_the_override_is_logged_with_BOTH_numbers(self):
        """A supplied APL far from the calibrated one is either the patch being wrong -- the case
        this exists for -- or the supplied value being wrong. The reader needs both to tell."""
        with self.assertLogs('pestifer.tasks.make_membrane_system', level='WARNING') as cm:
            _ov({'leaflet_apl': {'lower': 43.47}}, 41.32, 45.63)
        msg = '\n'.join(cm.output)
        self.assertIn('43.47', msg, 'the supplied value is not in the log')
        self.assertIn('45.63', msg, 'the DISPLACED calibrated value is not in the log')
        self.assertIn('-4.7%', msg, 'the log does not say how far off the patch was')
        self.assertNotIn('41.32', msg, 'the untouched leaflet should not be reported as overridden')


class TestTheCalibrationIsNotBypassed(unittest.TestCase):
    """NEGATIVE CONTROL for the thing that would make this feature harmful.

    If the override moved above the calibration, or short-circuited it, an unconverged patch
    would stop being reported and the override would be hiding exactly the defect it is usually
    reached for.
    """

    def setUp(self):
        self.src = inspect.getsource(MakeMembraneSystemTask.build_grid_membrane_asymmetric)

    def test_the_override_runs_AFTER_the_calibration_and_its_drift_check(self):
        i_cal = self.src.index('apl_upper = self.patchA.area')
        i_drift = self.src.index('area_drift')
        i_ov = self.src.index('_apply_leaflet_apl_override')
        self.assertLess(i_cal, i_drift, 'the drift check no longer follows the calibration')
        self.assertLess(i_drift, i_ov,
                        'the override runs before the convergence check, so a supplied APL now '
                        'suppresses the warning about an unequilibrated calibration cell')

    def test_the_calibration_is_still_computed_unconditionally(self):
        """Not guarded behind "unless an APL was supplied"."""
        for line in self.src.splitlines():
            if 'apl_upper = self.patchA.area' in line:
                self.assertEqual(line.strip().startswith('apl_upper'), True)
                self.assertEqual(len(line) - len(line.lstrip()), 8,
                                 f'calibration is nested under a conditional: {line!r}')
                break
        else:
            self.fail('the patch calibration has gone')


class TestTheSchemaDeclaresIt(unittest.TestCase):
    """A wrong-but-valid `type:` silently disables a whole subtree -- see CLAUDE.md."""

    def test_leaflet_apl_is_declared_a_dict_with_both_leaflets(self):
        import yaml
        from importlib.resources import files as pkg_files
        spec = yaml.safe_load(open(str(pkg_files('pestifer.schema').joinpath('base.yaml'))))

        def find(node, name):
            if isinstance(node, dict):
                if node.get('name') == name:
                    return node
                for v in node.values():
                    r = find(v, name)
                    if r is not None:
                        return r
            elif isinstance(node, list):
                for v in node:
                    r = find(v, name)
                    if r is not None:
                        return r
            return None

        node = find(spec, 'leaflet_apl')
        self.assertIsNotNone(node, 'leaflet_apl is not in the schema')
        self.assertEqual(node['type'], 'dict',
                         "declared anything but 'dict' and ycleptic walks straight past it")
        self.assertIn('text', node, 'a key with no text: crashes the RTD config-ref build')
        self.assertEqual(sorted(a['name'] for a in node['attributes']), ['lower', 'upper'])
        for a in node['attributes']:
            self.assertEqual(a['type'], 'float')


if __name__ == '__main__':
    unittest.main()
