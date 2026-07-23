"""Unit tests for pestifer.util.density_convergence (the NAMD-free convergence machinery)."""
import os
import tempfile
import unittest

import numpy as np

from pestifer.util.density_convergence import (
    ConvergenceParams,
    DensityConvergenceMonitor,
    next_chunk_steps,
    total_mass_amu,
    volume_to_density,
    xst_cell_volumes,
    xst_max_shrink_rate,
)


def _write_xst(path, rows):
    """rows: list of (ts, a_x, b_y, c_z) -> orthorhombic .xst lines (off-diagonal cell terms 0)."""
    with open(path, 'w') as f:
        f.write('# NAMD extended system trajectory\n')
        f.write('#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n')
        for ts, ax, by, cz in rows:
            f.write(f'{ts} {ax} 0 0  0 {by} 0  0 0 {cz}  0 0 0\n')


class TestXstParsing(unittest.TestCase):
    def test_volume_triple_product(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 't.xst')
            _write_xst(p, [(0, 10.0, 20.0, 30.0), (100, 9.0, 18.0, 27.0)])
            ts, vol = xst_cell_volumes(p)
            np.testing.assert_allclose(ts, [0, 100])
            np.testing.assert_allclose(vol, [10 * 20 * 30, 9 * 18 * 27])

    def test_empty_xst(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 't.xst')
            _write_xst(p, [])
            ts, vol = xst_cell_volumes(p)
            self.assertEqual(ts.size, 0)
            self.assertEqual(vol.size, 0)

    def test_density_conversion(self):
        # 1 amu in 1 A^3 -> 1.6605 g/cc
        d = volume_to_density(1.0, 1.0)
        self.assertAlmostEqual(float(d), 1.66053906660, places=6)
        # water-like: ~18 amu per ~30 A^3 -> ~1 g/cc
        d2 = float(volume_to_density(30.0, 18.0))
        self.assertTrue(0.9 < d2 < 1.1)


class TestShrinkRateAndChunk(unittest.TestCase):
    def test_shrink_rate_uses_most_shrunk_dim(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 't.xst')
            # over 100 steps: x shrinks 2 A, y shrinks 1, z grows; rate = 2/100
            _write_xst(p, [(0, 50.0, 50.0, 50.0), (100, 48.0, 49.0, 51.0)])
            self.assertAlmostEqual(xst_max_shrink_rate(p), 2.0 / 100)

    def test_shrink_rate_zero_when_growing(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 't.xst')
            _write_xst(p, [(0, 50.0, 50.0, 50.0), (100, 51.0, 51.0, 51.0)])
            self.assertEqual(xst_max_shrink_rate(p), 0.0)

    def test_shrink_rate_single_frame(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 't.xst')
            _write_xst(p, [(0, 50.0, 50.0, 50.0)])
            self.assertEqual(xst_max_shrink_rate(p), 0.0)

    def test_next_chunk_clamps(self):
        # settled box -> chunk_max
        self.assertEqual(next_chunk_steps(0.0, 4, 0.5, 500, 5000), 5000)
        # fast shrink -> throttled toward chunk_min
        # allowed = 0.5*4 / 0.1 = 20 steps -> clamped up to chunk_min
        self.assertEqual(next_chunk_steps(0.1, 4, 0.5, 500, 5000), 500)
        # moderate: allowed = 0.5*4/0.001 = 2000
        self.assertEqual(next_chunk_steps(0.001, 4, 0.5, 500, 5000), 2000)


class TestConvergenceMonitor(unittest.TestCase):
    def _params(self, **kw):
        base = dict(drift_tol=1e-3, precision_p=3.0, n_blocks=6, burn_in=100,
                    min_steps=200, n_consecutive=3)
        base.update(kw)
        return ConvergenceParams(**base)

    def test_below_min_steps_never_converges(self):
        mon = DensityConvergenceMonitor(self._params(min_steps=100000))
        t = np.arange(100, 5000, 100)
        mon.add_samples(t, np.full(t.size, 1.0))  # perfectly flat
        rep = mon.check()
        self.assertFalse(rep.converged)
        self.assertIn('min_steps', rep.reason)

    def test_flat_series_converges_after_hysteresis(self):
        p = self._params(min_steps=200, n_consecutive=3)
        mon = DensityConvergenceMonitor(p)
        # tiny deterministic ripple so std>0 but SEM/mean well below the gate
        step = 0
        results = []
        for _ in range(5):
            t = np.arange(step + 100, step + 2100, 100)
            d = 1.0 + 1e-6 * np.sin(t)  # essentially flat
            mon.add_samples(t, d)
            step += 2000
            results.append(mon.check())
        self.assertTrue(results[-1].converged)
        # must have taken at least n_consecutive passing checks
        self.assertGreaterEqual(results[-1].passes, 3)
        self.assertTrue(results[-1].precision_met)

    def test_drifting_series_does_not_converge(self):
        p = self._params(min_steps=200, n_consecutive=3)
        mon = DensityConvergenceMonitor(p)
        step = 0
        rep = None
        for _ in range(6):
            t = np.arange(step + 100, step + 2100, 100)
            d = 1.0 + 5e-4 * (t / 1000.0)  # steady ~0.05%/1000-step upward drift
            mon.add_samples(t, d)
            step += 2000
            rep = mon.check()
        self.assertFalse(rep.converged)
        self.assertGreater(rep.drift, p.drift_tol)

    def test_nan_triggers_blowup(self):
        mon = DensityConvergenceMonitor(self._params())
        t = np.arange(100, 3000, 100)
        d = np.full(t.size, 1.0)
        d[-1] = np.nan
        mon.add_samples(t, d)
        rep = mon.check()
        self.assertTrue(rep.blowup)
        self.assertFalse(rep.converged)

    def test_hysteresis_resets_on_failure(self):
        # a passing check followed by a fail must reset the consecutive counter
        p = self._params(min_steps=200, n_consecutive=2)
        mon = DensityConvergenceMonitor(p)
        t1 = np.arange(100, 2100, 100)
        mon.add_samples(t1, 1.0 + 1e-6 * np.sin(t1))
        r1 = mon.check()
        self.assertEqual(r1.passes, 1)
        # inject a drifting chunk -> should fail and reset
        t2 = np.arange(2100, 4100, 100)
        mon.add_samples(t2, 1.0 + 5e-3 * (t2 / 1000.0))
        r2 = mon.check()
        self.assertEqual(r2.passes, 0)
        self.assertFalse(r2.converged)


class TestTotalMass(unittest.TestCase):
    def test_total_mass_amu(self):
        psf = ('PSF\n\n       2 !NATOM\n'
               '       1 A    1    WAT  OH2  OT    -0.834000       15.9994           0\n'
               '       2 A    1    WAT  H1   HT     0.417000        1.0080           0\n')
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 't.psf')
            with open(p, 'w') as f:
                f.write(psf)
            self.assertAlmostEqual(total_mass_amu(p), 15.9994 + 1.008, places=4)


if __name__ == '__main__':
    unittest.main()
