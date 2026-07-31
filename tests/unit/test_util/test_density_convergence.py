"""Unit tests for pestifer.util.density_convergence (the NAMD-free convergence machinery)."""
import os
import tempfile
import unittest

import numpy as np

from pathlib import Path

from pestifer.util.density_convergence import (
    ConvergenceParams,
    DensityConvergenceMonitor,
    JointConvergence,
    is_patch_grid_crash,
    membrane_leaflet_geometry,
    next_chunk_steps,
    integrated_autocorr_time,
    parse_patch_grid,
    quantize_steps,
    total_atoms,
    total_mass_amu,
    volume_to_density,
    xst_cell_areas,
    xst_cell_volumes,
    xst_max_shrink_rate,
)
from pestifer.util.density_convergence import _read_coor_xyz


def _write_xst(path, rows):
    """rows: list of (ts, a_x, b_y, c_z) -> orthorhombic .xst lines (off-diagonal cell terms 0)."""
    with open(path, 'w') as f:
        f.write('# NAMD extended system trajectory\n')
        f.write('#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n')
        for ts, ax, by, cz in rows:
            f.write(f'{ts} {ax} 0 0  0 {by} 0  0 0 {cz}  0 0 0\n')


def _write_psf_pdb(psf_path, pdb_path, atoms):
    """atoms: list of (segname, resid, resname, name, mass, x, y, z) -> a minimal matched PSF+PDB."""
    with open(psf_path, 'w') as f:
        f.write('PSF CMAP\n\n       1 !NTITLE\n REMARKS synthetic\n\n')
        f.write(f'{len(atoms):8d} !NATOM\n')
        for i, (seg, resid, resn, name, mass, x, y, z) in enumerate(atoms, start=1):
            f.write(f'{i:8d} {seg:<4s} {resid:<4s} {resn:<4s} {name:<4s} {name:<4s} '
                    f'0.000000 {mass:10.4f}           0\n')
    with open(pdb_path, 'w') as f:
        for i, (seg, resid, resn, name, mass, x, y, z) in enumerate(atoms, start=1):
            f.write(f'ATOM  {i:5d} {name:<4s} {resn:<4s}{" ":1s}{resid:>4s}    '
                    f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n')


class TestXstParsing(unittest.TestCase):
    def test_volume_triple_product(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 't.xst')
            _write_xst(p, [(0, 10.0, 20.0, 30.0), (100, 9.0, 18.0, 27.0)])
            ts, vol = xst_cell_volumes(p)
            np.testing.assert_allclose(ts, [0, 100])
            np.testing.assert_allclose(vol, [10 * 20 * 30, 9 * 18 * 27])

    def test_lateral_area_cross_product(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 't.xst')
            _write_xst(p, [(0, 10.0, 20.0, 30.0), (100, 9.0, 18.0, 27.0)])
            ts, area = xst_cell_areas(p)
            np.testing.assert_allclose(ts, [0, 100])
            np.testing.assert_allclose(area, [10 * 20, 9 * 18])  # |a x b| = a_x*b_y (orthorhombic)

    def test_empty_xst(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 't.xst')
            _write_xst(p, [])
            ts, vol = xst_cell_volumes(p)
            self.assertEqual(ts.size, 0)
            self.assertEqual(vol.size, 0)
            self.assertEqual(xst_cell_areas(p)[1].size, 0)

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

    def test_quantize_steps(self):
        # NAMD `run` counts must be whole cycles: round DOWN to a multiple of steps_per_cycle
        self.assertEqual(quantize_steps(675, 10), 670)
        self.assertEqual(quantize_steps(200, 10), 200)
        self.assertEqual(quantize_steps(7, 10), 10)        # floors at one cycle
        self.assertEqual(quantize_steps(675, 50), 650)
        # minimum is rounded UP to a whole number of cycles, then used as the floor
        self.assertEqual(quantize_steps(5, 10, minimum=20), 20)
        self.assertEqual(quantize_steps(5, 10, minimum=15), 20)
        self.assertEqual(quantize_steps(1000, 10, minimum=20), 1000)

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
        base = dict(drift_tol=1e-3, precision_p=3.0, burn_in=100,
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

    def test_signed_drift_carries_direction(self):
        # signed_drift keeps the trend sign (>0 growing, <0 shrinking); drift is its magnitude
        p = self._params(min_steps=200, n_consecutive=3)
        for sign in (+1.0, -1.0):
            mon = DensityConvergenceMonitor(p)
            step, rep = 0, None
            for _ in range(4):
                t = np.arange(step + 100, step + 2100, 100)
                d = 1.0 + sign * 5e-4 * (t / 1000.0)   # steady up (or down) drift
                mon.add_samples(t, d)
                step += 2000
                rep = mon.check()
            self.assertEqual(rep.signed_drift > 0, sign > 0)   # direction preserved
            self.assertAlmostEqual(rep.drift, abs(rep.signed_drift))   # magnitude drives the gate

    def test_trailing_window_ages_out_transient(self):
        # A ramp-then-plateau series: a full-history window (window_frac=1) keeps seeing the ramp and
        # never converges; a trailing window (0.5) lets it age out and converges once the trailing
        # half is all plateau.  Mirrors the BPTI validation finding.
        rng_t = np.arange(100, 40100, 100)
        # ramp 0.95 -> 1.03 over the first 10000 steps, then flat 1.03 with tiny ripple
        plateau = 1.03
        d = np.where(rng_t < 10000,
                     0.95 + (plateau - 0.95) * (rng_t / 10000.0),
                     plateau + 1e-5 * np.sin(rng_t))

        def run(window_frac):
            mon = DensityConvergenceMonitor(self._params(min_steps=200, n_consecutive=2,
                                                         burn_in=1000))
            mon.params.window_frac = window_frac
            converged_at = None
            # feed in 2000-step chunks, checking at each boundary
            for lo in range(0, rng_t.size, 20):
                sl = slice(lo, lo + 20)
                if not len(rng_t[sl]):
                    break
                mon.add_samples(rng_t[sl], d[sl])
                r = mon.check()
                if r.converged and converged_at is None:
                    converged_at = r.last_step
            return converged_at

        self.assertIsNone(run(1.0))          # full history: transient never leaves -> never converges
        self.assertIsNotNone(run(0.5))       # trailing half: converges once past the ramp
        self.assertGreater(run(0.5), 15000)  # ...but only after the trailing window clears the 10k ramp

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
        p.window_frac = 1.0
        mon = DensityConvergenceMonitor(p)
        t1 = np.arange(300, 300 + 100 * 60, 100)  # 60 flat samples -> a legitimate reliable pass
        mon.add_samples(t1, 1.0 + 1e-6 * np.sin(t1))
        r1 = mon.check()
        self.assertEqual(r1.passes, 1)
        # inject a strongly drifting chunk -> should fail on drift and reset
        t2 = np.arange(t1[-1] + 100, t1[-1] + 100 + 100 * 30, 100)
        mon.add_samples(t2, 1.0 + 5e-3 * (t2 / 1000.0))
        r2 = mon.check()
        self.assertEqual(r2.passes, 0)
        self.assertFalse(r2.converged)

    def test_burn_in_is_relative_to_first_sample(self):
        # The burn-in must be measured from the first sample, not as an absolute timestep: a run that
        # continues from a high `firsttimestep` (e.g. 8000, after NVT warm-up) must still discard its
        # own leading transient.  Two identical ramp-then-plateau series that differ only by an absolute
        # time offset must therefore give the same verdict at every step.
        def verdicts(offset):
            p = self._params(min_steps=200, n_consecutive=2, burn_in=1000)
            p.window_frac = 0.5
            mon = DensityConvergenceMonitor(p)
            rel = np.arange(100, 20100, 100)
            d = np.where(rel < 5000, 0.95 + (1.03 - 0.95) * (rel / 5000.0),
                         1.03 + 1e-5 * np.sin(rel))
            out = []
            for lo in range(0, rel.size, 20):
                sl = slice(lo, lo + 20)
                if not len(rel[sl]):
                    break
                mon.add_samples(rel[sl] + offset, d[sl])
                out.append(mon.check().converged)
            return out
        # absolute offset of 8000 must not change the convergence trajectory
        self.assertEqual(verdicts(0), verdicts(8000))
        self.assertTrue(any(verdicts(8000)))  # and it does converge

    def test_autocorr_corrected_sem_exceeds_naive(self):
        # For a correlated (but flat) series the honest SEM (sigma/sqrt(N_eff), N_eff=N/tau) must be
        # LARGER than the naive sigma/sqrt(N), because tau>1 shrinks the effective sample count.
        rng = np.random.default_rng(7)
        # AR(1) with phi=0.8 -> tau_int = (1+phi)/(1-phi) = 9
        n = 4000
        x = np.zeros(n)
        for i in range(1, n):
            x[i] = 0.8 * x[i - 1] + rng.standard_normal()
        tau = integrated_autocorr_time(x)
        self.assertGreater(tau, 4.0)        # clearly correlated (true tau_int = 9)
        sigma = x.std(ddof=1)
        naive_sem = sigma / np.sqrt(n)
        honest_sem = sigma / np.sqrt(n / tau)
        self.assertGreater(honest_sem, 2 * naive_sem)   # ~3x for tau=9

    def test_integrated_autocorr_time_white_vs_correlated(self):
        rng = np.random.default_rng(11)
        white = rng.standard_normal(5000)
        self.assertLess(integrated_autocorr_time(white), 2.0)   # ~1 for independent samples
        # smooth (integrated white noise, then differenced lightly) -> strongly correlated
        phi = 0.9
        x = np.zeros(5000)
        for i in range(1, x.size):
            x[i] = phi * x[i - 1] + rng.standard_normal()
        self.assertGreater(integrated_autocorr_time(x), 8.0)    # true tau_int = 19
        self.assertEqual(integrated_autocorr_time(np.full(50, 3.0)), 1.0)  # zero variance -> 1

    def test_short_window_not_trusted_until_several_taus(self):
        # A correlated series must NOT be declared converged while the window spans too few
        # autocorrelation times, even if it happens to look flat -- tau is not estimable there.
        p = self._params(min_steps=200, n_consecutive=1, autocorr_reliability=6.0)
        p.window_frac = 1.0
        mon = DensityConvergenceMonitor(p)
        rng = np.random.default_rng(5)
        # strongly correlated flat series
        x = np.zeros(400)
        for i in range(1, x.size):
            x[i] = 0.85 * x[i - 1] + rng.standard_normal()
        t = np.arange(100, 100 + 100 * x.size, 100)
        d = 1.0 + 1e-4 * x
        # feed the first few samples only: window too short vs tau -> must not converge / not reliable
        mon.add_samples(t[:20], d[:20])
        r = mon.check()
        self.assertFalse(r.converged)
        self.assertIn('autocorrelation', r.reason)
        self.assertIsNotNone(r.tau_int)

    def test_drift_uses_upper_confidence_bound(self):
        # Convergence gates on the UPPER confidence bound of |drift|, not the point estimate.  On the
        # same flat data, a modest drift_conf certifies while an absurd one refuses -- proving the bound
        # (not the raw slope) is what's tested.  And drift_hi >= drift always.
        t = np.arange(300, 300 + 100 * 40, 100)
        rng = np.random.default_rng(9)
        d = 1.0 + 5e-4 * rng.standard_normal(t.size)   # flat, small noise
        p1 = self._params(min_steps=200, n_consecutive=1, drift_tol=3e-3, drift_conf=2.0)
        p1.window_frac = 1.0
        m1 = DensityConvergenceMonitor(p1); m1.add_samples(t, d); r1 = m1.check()
        self.assertGreaterEqual(r1.drift_hi, r1.drift)
        self.assertTrue(r1.converged)
        p2 = self._params(min_steps=200, n_consecutive=1, drift_tol=3e-3, drift_conf=1000.0)
        p2.window_frac = 1.0
        m2 = DensityConvergenceMonitor(p2); m2.add_samples(t, d); r2 = m2.check()
        self.assertFalse(r2.converged)
        self.assertIn('confidently', r2.reason)


class TestPatchGridCrash(unittest.TestCase):
    _CRASH = ('Info: Entering startup phase 0\n'
              'FATAL ERROR: Periodic cell has become too small for original patch grid!\n'
              'Possible solutions are to restart from a recent checkpoint,\n'
              'increase margin, or disable useFlexibleCell for liquid simulation.\n')

    def test_detects_patch_grid_crash(self):
        self.assertTrue(is_patch_grid_crash(self._CRASH))

    def test_ignores_other_failures(self):
        self.assertFalse(is_patch_grid_crash('FATAL ERROR: Constraint failure in RATTLE algorithm!'))
        self.assertFalse(is_patch_grid_crash(''))
        self.assertFalse(is_patch_grid_crash(None))

    def test_case_insensitive(self):
        self.assertTrue(is_patch_grid_crash('periodic cell too small for the PATCH GRID'))


class TestParsePatchGrid(unittest.TestCase):
    _LOG = (
        'Info: CUTOFF                10\n'
        'Info: PAIRLIST DISTANCE     11.5\n'
        'Info: MARGIN                4\n'
        'Info: PERIODIC CELL BASIS 1  62 0 0\n'
        'Info: PERIODIC CELL BASIS 2  0 62 0\n'
        'Info: PERIODIC CELL BASIS 3  0 0 124\n'
        'Info: PATCH GRID IS 5 (PERIODIC) BY 5 (PERIODIC) BY 10 (PERIODIC)\n'
    )

    def test_parses_grid_and_derived(self):
        info = parse_patch_grid(self._LOG)
        self.assertEqual(info['patches'], (5, 5, 10))
        self.assertEqual(info['cell_edges'], (62.0, 62.0, 124.0))
        self.assertEqual(info['cutoff'], 10.0)
        self.assertEqual(info['pairlistdist'], 11.5)
        self.assertEqual(info['margin'], 4.0)
        # min patch dim = 62/5 = 12.4 (both x and z give 12.4); headroom to pairlist = 12.4-11.5=0.9
        self.assertAlmostEqual(info['min_patch_dim'], 12.4, places=3)
        self.assertAlmostEqual(info['shrink_headroom'], 0.9, places=3)

    def test_no_grid_returns_none(self):
        self.assertIsNone(parse_patch_grid('Info: nothing useful here'))
        self.assertIsNone(parse_patch_grid(''))

    def test_grid_without_cell(self):
        info = parse_patch_grid('Info: PATCH GRID IS 3 BY 4 BY 5\n')
        self.assertEqual(info['patches'], (3, 4, 5))
        self.assertNotIn('min_patch_dim', info)

    def test_gpu_few_patch_omits_negative_headroom(self):
        # GPU-resident NAMD tiles a small box into few patches, so a patch dimension can be BELOW the
        # pairlist distance while NAMD runs fine.  min_patch_dim is still reported (a real fact), but
        # the derived "headroom" must be omitted rather than reported negative/misleading.
        gpu_log = (
            'Info: CUTOFF                10\n'
            'Info: PAIRLIST DISTANCE     11.5\n'
            'Info: PERIODIC CELL BASIS 1  58 0 0\n'
            'Info: PERIODIC CELL BASIS 2  0 51 0\n'
            'Info: PERIODIC CELL BASIS 3  0 0 52\n'
            'Info: PATCH GRID IS 6 (PERIODIC) BY 2 (PERIODIC) BY 2 (PERIODIC)\n'
        )
        info = parse_patch_grid(gpu_log)
        self.assertEqual(info['patches'], (6, 2, 2))
        # x: 58/6 = 9.67 < pairlistdist 11.5 -> non-positive headroom
        self.assertAlmostEqual(info['min_patch_dim'], 58 / 6, places=3)
        self.assertNotIn('shrink_headroom', info)


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
            self.assertEqual(total_atoms(p), 2)


class TestJointConvergence(unittest.TestCase):
    def _mon(self, **kw):
        base = dict(drift_tol=1e-3, precision_p=3.0, burn_in=100, min_steps=200)
        base.update(kw)
        m = DensityConvergenceMonitor(ConvergenceParams(**base))
        m.params.window_frac = 1.0
        return m

    def _flat(self, mon, seed):
        rng = np.random.default_rng(seed)
        for lo in range(80):
            t = np.arange(lo * 2000 + 100, lo * 2000 + 2100, 100)
            mon.add_samples(t, 1.0 + 1e-5 * rng.standard_normal(t.size))

    def test_both_must_converge_jointly(self):
        # density flat (would converge alone) but area still ramping -> joint does NOT converge.
        d = self._mon(); a = self._mon()
        jc = JointConvergence({'density': d, 'area': a}, n_consecutive=2)
        rng = np.random.default_rng(0)
        conv = False
        for lo in range(60):
            t = np.arange(lo * 2000 + 100, lo * 2000 + 2100, 100)
            jc.add_samples('density', t, 1.0 + 1e-5 * rng.standard_normal(t.size))       # flat
            jc.add_samples('area', t, 100.0 + 5e-3 * t)                                   # steady ramp
            if jc.check().converged:
                conv = True; break
        self.assertFalse(conv)
        # the joint report names the lagging observable
        self.assertIn('area', jc.check().reason)

    def test_converges_when_all_flat(self):
        d = self._mon(); a = self._mon()
        jc = JointConvergence({'density': d, 'area': a}, n_consecutive=2)
        rd = np.random.default_rng(1); ra = np.random.default_rng(2)
        conv_step = None
        for lo in range(80):
            t = np.arange(lo * 2000 + 100, lo * 2000 + 2100, 100)
            jc.add_samples('density', t, 1.0 + 1e-5 * rd.standard_normal(t.size))
            jc.add_samples('area', t, 100.0 + 1e-3 * ra.standard_normal(t.size))
            r = jc.check()
            if r.converged:
                conv_step = r.last_step; break
        self.assertIsNotNone(conv_step)
        self.assertGreaterEqual(jc._passes, 2)

    def test_blowup_in_one_observable(self):
        d = self._mon(); a = self._mon()
        jc = JointConvergence({'density': d, 'area': a}, n_consecutive=2)
        t = np.arange(100, 3000, 100)
        jc.add_samples('density', t, np.full(t.size, 1.0))
        av = np.full(t.size, 100.0); av[-1] = np.inf
        jc.add_samples('area', t, av)
        r = jc.check()
        self.assertTrue(r.blowup)
        self.assertFalse(r.converged)
        self.assertIn('area', r.reason)

    def test_sub_monitors_forced_single_check_hysteresis(self):
        # JointConvergence owns the hysteresis: it sets each monitor's n_consecutive to 1.
        d = self._mon(n_consecutive=3)
        JointConvergence({'density': d}, n_consecutive=4)
        self.assertEqual(d.params.n_consecutive, 1)


if __name__ == '__main__':
    unittest.main()


class TestMembraneLeafletGeometry(unittest.TestCase):
    """Per-leaflet lipid counts + protein-corrected APL from an embedded membrane frame."""

    def _synthetic(self):
        # midplane at z=0; 3 lipids per leaflet (mass-weighted mean-z puts each on its side);
        # protein: a 10x10 square in the upper leaflet (hull area 100) and a 5x5 in the lower (25)
        atoms = []
        for r, (px, py) in enumerate([(0, 0), (30, 0), (0, 30)], start=1):
            atoms.append(('MEMB', str(r), 'DMPC', 'P', 31.0, px, py, -20.0))
            atoms.append(('MEMB', str(r), 'DMPC', 'C22', 12.0, px, py, -5.0))
        for r, (px, py) in enumerate([(0, 0), (30, 0), (0, 30)], start=4):
            atoms.append(('MEMB', str(r), 'DMPC', 'P', 31.0, px, py, 20.0))
            atoms.append(('MEMB', str(r), 'DMPC', 'C22', 12.0, px, py, 5.0))
        for (x, y) in [(0, 0), (10, 0), (10, 10), (0, 10)]:      # upper protein square -> 100
            atoms.append(('PROT', '1', 'ALA', 'CA', 12.0, x, y, 10.0))
        for (x, y) in [(0, 0), (5, 0), (5, 5), (0, 5)]:          # lower protein square -> 25
            atoms.append(('PROT', '2', 'ALA', 'CA', 12.0, x, y, -10.0))
        atoms.append(('WT1', '1', 'TIP3', 'OH2', 16.0, 50.0, 50.0, 0.0))   # ignored
        return atoms

    def test_counts_footprints_and_apl(self):
        atoms = self._synthetic()
        with tempfile.TemporaryDirectory() as d:
            psf, pdb = os.path.join(d, 's.psf'), os.path.join(d, 's.pdb')
            _write_psf_pdb(psf, pdb, atoms)
            g = membrane_leaflet_geometry(psf, pdb)
        self.assertEqual((g.n_lower, g.n_upper), (3, 3))
        self.assertEqual(g.n_lipids_total, 6)
        self.assertAlmostEqual(g.midplane_z, 0.0, places=6)
        self.assertAlmostEqual(g.a_prot_upper, 100.0, places=4)
        self.assertAlmostEqual(g.a_prot_lower, 25.0, places=4)
        # protein-corrected APL is below the naive area/n_lipids, by footprint/n_lipids
        A = 600.0
        self.assertAlmostEqual(g.apl_upper(A), (600.0 - 100.0) / 3, places=4)
        self.assertAlmostEqual(g.apl_lower(A), (600.0 - 25.0) / 3, places=4)
        self.assertLess(g.apl_mean(A), A / 3)   # naive per-leaflet would be 200

    def test_degenerate_footprint_is_zero(self):
        # drop the lower protein square -> fewer than 3 atoms in that leaflet -> 0 footprint
        atoms = [a for a in self._synthetic() if not (a[0] == 'PROT' and a[1] == '2')]
        with tempfile.TemporaryDirectory() as d:
            psf, pdb = os.path.join(d, 's.psf'), os.path.join(d, 's.pdb')
            _write_psf_pdb(psf, pdb, atoms)
            g = membrane_leaflet_geometry(psf, pdb)
        self.assertEqual(g.a_prot_lower, 0.0)
        self.assertAlmostEqual(g.apl_lower(600.0), 200.0, places=4)   # falls back to naive when no protein

    def test_no_lipids_raises(self):
        atoms = [('PROT', '1', 'ALA', 'CA', 12.0, 0.0, 0.0, 0.0),
                 ('WT1', '1', 'TIP3', 'OH2', 16.0, 5.0, 5.0, 5.0)]
        with tempfile.TemporaryDirectory() as d:
            psf, pdb = os.path.join(d, 's.psf'), os.path.join(d, 's.pdb')
            _write_psf_pdb(psf, pdb, atoms)
            with self.assertRaises(ValueError):
                membrane_leaflet_geometry(psf, pdb)

    def test_accepts_pathlib_coord(self):
        # regression: coord artifacts are pathlib.Path; _read_coor_xyz used to call str-only .lower()
        # on the path and blow up, silently forcing the lipids_per_leaflet fallback
        atoms = self._synthetic()
        with tempfile.TemporaryDirectory() as d:
            psf, pdb = Path(d) / 's.psf', Path(d) / 's.pdb'
            _write_psf_pdb(str(psf), str(pdb), atoms)
            xyz = _read_coor_xyz(pdb, len(atoms))             # Path in, no AttributeError
            self.assertEqual(xyz.shape[1], 3)
            g = membrane_leaflet_geometry(psf, pdb)           # end-to-end with Path args
        self.assertEqual((g.n_lower, g.n_upper), (3, 3))


class TestMonitorReset(unittest.TestCase):
    def test_monitor_reset_clears_samples_and_passes(self):
        mon = DensityConvergenceMonitor(ConvergenceParams(min_steps=0, n_consecutive=1))
        t = np.arange(100, 100 + 100 * 40, 100)
        mon.add_samples(t, 1.0 + 1e-6 * np.sin(t))
        mon.check()
        self.assertTrue(mon._t)
        mon.reset()
        self.assertEqual(mon._t, [])
        self.assertEqual(mon._d, [])
        self.assertEqual(mon._passes, 0)

    def test_joint_reset_clears_all_monitors(self):
        a = DensityConvergenceMonitor(ConvergenceParams(min_steps=0, n_consecutive=1))
        b = DensityConvergenceMonitor(ConvergenceParams(min_steps=0, n_consecutive=1))
        jc = JointConvergence({'a': a, 'b': b}, n_consecutive=1)
        jc.add_samples('a', [100, 200], [1.0, 1.0])
        jc.add_samples('b', [100, 200], [1.0, 1.0])
        jc._passes = 2
        jc.reset()
        self.assertEqual(a._t, [])
        self.assertEqual(b._t, [])
        self.assertEqual(jc._passes, 0)
