"""Unit tests for MembraneEquilibrateTask's observable-parameter resolution and APL (no pipeline)."""
import unittest

from pestifer.tasks.membrane_equilibrate import MembraneEquilibrateTask
from pestifer.tasks.equilibrate_base import ChunkedEquilibrateTask
from pestifer.util.density_convergence import LeafletGeometry


def _task(specs, geom=None):
    # build a bare instance without running MDTask.__init__ -- _params_for reads self.specs, and the
    # APL helpers read _lipids_per_leaflet + _geom resolved in _setup (set them directly here).
    t = MembraneEquilibrateTask.__new__(MembraneEquilibrateTask)
    t.specs = specs
    t._lipids_per_leaflet = int(specs.get('lipids_per_leaflet') or 0)
    t._geom = geom
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

    def test_apl_spec_fallback(self):
        # no measured geometry: fall back to naive area/lipids_per_leaflet (same for both leaflets)
        self.assertEqual(_task(dict(BASE, lipids_per_leaflet=0))._apl_pair(4200.0), (None, None))
        self.assertIsNone(_task(dict(BASE, lipids_per_leaflet=0))._apl_mean(4200.0))
        lo, up = _task(dict(BASE, lipids_per_leaflet=70))._apl_pair(4200.0)
        self.assertAlmostEqual(lo, 60.0)
        self.assertAlmostEqual(up, 60.0)
        self.assertAlmostEqual(_task(dict(BASE, lipids_per_leaflet=70))._apl_mean(4200.0), 60.0)
        self.assertEqual(_task(dict(BASE, lipids_per_leaflet=70))._apl_pair(None), (None, None))

    def test_apl_protein_corrected_per_leaflet(self):
        # measured geometry: APL is (area - protein_footprint)/n_lipids, per leaflet
        geom = LeafletGeometry(n_lower=70, n_upper=70, a_prot_lower=140.0, a_prot_upper=700.0,
                               midplane_z=0.0, n_lipids_total=140)
        t = _task(dict(BASE, lipids_per_leaflet=0), geom=geom)   # geometry wins over (absent) spec
        lo, up = t._apl_pair(4200.0)
        self.assertAlmostEqual(lo, (4200.0 - 140.0) / 70)        # 58.0
        self.assertAlmostEqual(up, (4200.0 - 700.0) / 70)        # 50.0
        self.assertAlmostEqual(t._apl_mean(4200.0), 54.0)
        # correction lowers APL vs the naive 60.0
        self.assertLess(t._apl_mean(4200.0), 60.0)


class TestSubclassing(unittest.TestCase):
    def test_is_chunked_equilibrate(self):
        self.assertTrue(issubclass(MembraneEquilibrateTask, ChunkedEquilibrateTask))
        self.assertEqual(MembraneEquilibrateTask._ensemble, 'NPgT')


if __name__ == '__main__':
    unittest.main()
