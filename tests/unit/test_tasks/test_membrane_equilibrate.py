"""Unit tests for MembraneEquilibrateTask's observable-parameter resolution and APL (no pipeline)."""
import unittest

from pestifer.tasks.membrane_equilibrate import MembraneEquilibrateTask
from pestifer.tasks.equilibrate_base import ChunkedEquilibrateTask


def _task(specs):
    # build a bare instance without running MDTask.__init__ -- _params_for reads self.specs, and _apl
    # reads the _lipids_per_leaflet resolved in _setup (set it directly here).
    t = MembraneEquilibrateTask.__new__(MembraneEquilibrateTask)
    t.specs = specs
    t._lipids_per_leaflet = int(specs.get('lipids_per_leaflet') or 0)
    return t


BASE = dict(drift_tol=0.002, precision_p=3.0, drift_conf=0.0, burn_in=2000,
            window_frac=0.5, n_consecutive=3, autocorr_reliability=6.0)


class TestParamResolution(unittest.TestCase):
    def test_area_defaults_to_density(self):
        t = _task(dict(BASE))
        dp = t._params_for('density', 4000)
        ap = t._params_for('area', 4000)
        # with no area_* overrides, area mirrors density
        self.assertEqual(ap.drift_tol, dp.drift_tol)
        self.assertEqual(ap.precision_p, dp.precision_p)
        self.assertEqual(ap.burn_in, dp.burn_in)
        self.assertEqual(ap.min_steps, 4000)

    def test_area_override_applies(self):
        specs = dict(BASE, area_drift_tol=0.004, area_precision_p=2.0)
        t = _task(specs)
        dp = t._params_for('density', 4000)
        ap = t._params_for('area', 4000)
        self.assertEqual(dp.drift_tol, 0.002)          # density unchanged
        self.assertEqual(ap.drift_tol, 0.004)          # area overridden
        self.assertEqual(ap.precision_p, 2.0)
        self.assertEqual(ap.burn_in, dp.burn_in)       # unspecified area key still mirrors density

    def test_apl(self):
        self.assertIsNone(_task(dict(BASE, lipids_per_leaflet=0))._apl(4200.0))
        self.assertAlmostEqual(_task(dict(BASE, lipids_per_leaflet=70))._apl(4200.0), 60.0)
        self.assertIsNone(_task(dict(BASE, lipids_per_leaflet=70))._apl(None))


class TestSubclassing(unittest.TestCase):
    def test_is_chunked_equilibrate(self):
        self.assertTrue(issubclass(MembraneEquilibrateTask, ChunkedEquilibrateTask))
        self.assertEqual(MembraneEquilibrateTask._ensemble, 'NPgT')


if __name__ == '__main__':
    unittest.main()
