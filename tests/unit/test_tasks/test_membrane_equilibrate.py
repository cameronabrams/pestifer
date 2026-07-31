"""Unit tests for MembraneEquilibrateTask's observable-parameter resolution and APL (no pipeline)."""
import unittest
from unittest import mock

import numpy as np

from pestifer.tasks import membrane_equilibrate as me_mod
from pestifer.tasks.membrane_equilibrate import MembraneEquilibrateTask
from pestifer.tasks.equilibrate_base import ChunkedEquilibrateTask
from pestifer.util.density_convergence import (
    LeafletGeometry, DensityConvergenceMonitor, JointConvergence,
)


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


# --- two-stage barostat-mode + phase machine ------------------------------------------------------

# loose gates so a handful of flat chunks certifies convergence in-test (drift_tol huge, precision gate
# = drift_tol/precision_p huge, no autocorrelation-reliability floor, no burn-in, whole-history window)
LOOSE = dict(drift_tol=1.0, precision_p=0.01, drift_conf=0.0, burn_in=0,
             window_frac=1.0, autocorr_reliability=0.0, area_min_steps=0)


def _machine(two_stage=True, n_consecutive=1, user_other=None):
    """A MembraneEquilibrateTask with the post-_setup state built by hand (no file IO)."""
    t = MembraneEquilibrateTask.__new__(MembraneEquilibrateTask)
    t.taskname = 'membrane_equilibrate'
    t.specs = dict(LOOSE, n_consecutive=n_consecutive, two_stage=two_stage)
    if user_other is not None:
        t.specs['other_parameters'] = dict(user_other)
    t._mass_amu = 40000.0
    t._geom = None
    t._lipids_per_leaflet = 233
    t._min_steps = 0
    t._area_min = 0
    t._n_consecutive = n_consecutive
    t._mon_d = DensityConvergenceMonitor(t._params_for('density', 0))
    t._mon_a = DensityConvergenceMonitor(t._params_for('area', 0))
    t._all_t, t._all_d, t._all_a = [], [], []
    t._rows = []
    t._phase2_start_step = None
    t._user_other = dict(t.specs.get('other_parameters', {}))
    t._two_stage = two_stage
    t._phase = 1 if two_stage else 2
    t._set_barostat_mode(t._phase)
    if two_stage:
        t._conv = JointConvergence({'density': t._mon_d}, n_consecutive=n_consecutive)
    else:
        t._conv = JointConvergence({'density': t._mon_d, 'area': t._mon_a}, n_consecutive=n_consecutive)
    return t


def _flat_chunk(t, start_step, vol=1.0e6, area=1.18e4, n=30, step=100):
    """Feed one chunk of ~flat samples through _ingest_chunk with the .xst readers stubbed."""
    ts = np.arange(start_step, start_step + n * step, step, dtype=float)
    # tiny deterministic jitter so std>0 (keeps tau_int well-defined) but drift ~0
    jit = 1.0 + 1e-6 * np.sin(np.arange(n))
    vols = np.full(n, vol) * jit
    areas = np.full(n, area) * jit
    with mock.patch.object(me_mod, 'xst_cell_volumes', return_value=(ts, vols)), \
         mock.patch.object(me_mod, 'xst_cell_areas', return_value=(ts, areas)):
        return t._ingest_chunk('ignored.xst', int(ts[-1]), n * step, len(t._rows) + 1)


class TestBarostatMode(unittest.TestCase):
    def test_phase1_constant_area(self):
        t = _machine(two_stage=True)
        op = t.specs['other_parameters']
        self.assertEqual(op['useconstantarea'], 'yes')
        self.assertEqual(op['useconstantratio'], 'no')

    def test_phase2_tensionless(self):
        t = _machine(two_stage=False)   # single-stage starts in phase 2
        op = t.specs['other_parameters']
        self.assertEqual(op['useconstantarea'], 'no')
        self.assertEqual(op['useconstantratio'], 'yes')

    def test_user_other_parameters_preserved(self):
        t = _machine(two_stage=True, user_other={'watchdog': 'on'})
        op = t.specs['other_parameters']
        self.assertEqual(op['watchdog'], 'on')          # user key survives the merge
        self.assertEqual(op['useconstantarea'], 'yes')  # mode keys added on top


class TestPhaseMachine(unittest.TestCase):
    def test_stage1_gates_density_only_then_hands_off(self):
        t = _machine(two_stage=True, n_consecutive=1)
        blowup, converged = _flat_chunk(t, 0)
        # density settled at constant area -> hand off, but the TASK is not done yet
        self.assertFalse(converged)
        self.assertFalse(blowup)
        self.assertEqual(t._phase, 2)
        self.assertEqual(t._phase2_start_step, int(t._rows[-1][1]))
        # params flipped to tensionless, monitors reset for stage-2-only assessment
        self.assertEqual(t.specs['other_parameters']['useconstantratio'], 'yes')
        self.assertEqual(t.specs['other_parameters']['useconstantarea'], 'no')
        self.assertEqual(t._mon_d._t, [])
        self.assertEqual(t._mon_a._t, [])

    def test_stage2_gates_area_and_converges(self):
        t = _machine(two_stage=True, n_consecutive=1)
        _flat_chunk(t, 0)                       # stage 1 -> handoff (phase now 2)
        start = t._phase2_start_step
        blowup, converged = _flat_chunk(t, start)   # stage 2: area now a live observable
        self.assertFalse(blowup)
        self.assertTrue(converged)              # both density and area stationary -> done
        # stage 2 fed the area monitor (stage 1 did not)
        self.assertGreater(len(t._mon_a._t), 0)

    def test_stage1_does_not_feed_area_monitor(self):
        t = _machine(two_stage=True, n_consecutive=2)   # keep it in stage 1 (need 2 passes)
        _flat_chunk(t, 0)
        # after one chunk, not yet handed off; area monitor must be empty (area is pinned)
        self.assertEqual(t._phase, 1)
        self.assertEqual(t._mon_a._t, [])
        self.assertGreater(len(t._mon_d._t), 0)

    def test_single_stage_converges_without_handoff(self):
        t = _machine(two_stage=False, n_consecutive=1)
        blowup, converged = _flat_chunk(t, 0)
        self.assertTrue(converged)              # joint density+area from the start
        self.assertIsNone(t._phase2_start_step) # no stage boundary
        self.assertEqual(t._phase, 2)

    def test_rows_carry_stage(self):
        t = _machine(two_stage=True, n_consecutive=1)
        _flat_chunk(t, 0)
        self.assertEqual(t._rows[-1][4], 1)     # the handoff chunk is still tagged stage 1


if __name__ == '__main__':
    unittest.main()
