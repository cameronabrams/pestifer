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
    t._area_plateau_tol = float(t.specs.get('area_plateau_tol') or 0.0)
    t._last_plateau_drift = None
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


class TestPlateauGate:
    """The cumulative area-plateau gate (area_plateau_tol) for calibration convergence."""

    def _series(self, t, areas, start=0, step=100, phase2_start=0):
        import numpy as np
        n = len(areas)
        t._all_t = list(np.arange(start, start + n * step, step, dtype=float))
        t._all_a = list(map(float, areas))
        t._phase2_start_step = phase2_start

    def test_quarter_drift_flat_is_zero(self):
        import numpy as np
        t = _machine(two_stage=False)
        self._series(t, np.full(200, 4500.0))
        assert abs(t._area_plateau_drift()) < 1e-9

    def test_quarter_drift_detects_monotone_creep(self):
        import numpy as np
        t = _machine(two_stage=False)
        # steady descent 4600 -> 4500 over the series: 4th-quarter mean below 3rd-quarter mean
        self._series(t, np.linspace(4600.0, 4500.0, 200))
        pd = t._area_plateau_drift()
        assert pd is not None and pd < -0.005   # cumulative quarter-drift exceeds a 0.5% gate

    def test_quarter_drift_only_stage2(self):
        # stage-1 (pinned, flat) samples before phase2_start must be excluded
        import numpy as np
        t = _machine(two_stage=True)
        flat = np.full(100, 4700.0)              # stage 1: pinned
        creep = np.linspace(4600.0, 4500.0, 100) # stage 2: relaxing
        allA = np.concatenate([flat, creep])
        self._series(t, allA, phase2_start=100 * 100)   # phase2 starts after the 100 flat samples
        pd = t._area_plateau_drift()
        # computed over the stage-2 creep only (flat stage-1 excluded), so it sees the descent
        assert pd is not None and pd < -0.005

    def test_too_few_samples_returns_none(self):
        t = _machine(two_stage=False)
        self._series(t, [4500.0] * 6)            # tail < 8 samples
        assert t._area_plateau_drift() is None

    def test_gate_read_from_spec(self):
        import numpy as np
        t = _machine(two_stage=False)
        t.specs['area_plateau_tol'] = 0.005
        t._area_plateau_tol = 0.005
        # a flat area passes the gate; a creeping one would not (verified via _area_plateau_drift above)
        self._series(t, np.full(200, 4500.0))
        assert abs(t._area_plateau_drift()) < t._area_plateau_tol


# --- stop reasons and the convergence report -------------------------------------------------
#
# These are what a user reads and what run-record.json stores, so their content is part of the
# task's contract rather than decoration.  They were previously exercised only by running a real
# membrane equilibration.

import os
import tempfile

from pestifer.util.density_convergence import ConvergenceReport, JointConvergenceReport


def _report(mean=1.02, drift=1.0e-4, sem=1.0e-4, precision_met=True):
    return ConvergenceReport(mean_density=mean, signed_drift=drift, drift=abs(drift),
                             sem_over_mean=sem, precision_met=precision_met, passes=3)


def _joint(passes=3, converged=True, density=None, area=None, reason='ok'):
    return JointConvergenceReport(
        converged=converged, passes=passes,
        reports={'density': density if density is not None else _report(),
                 'area': area if area is not None else _report(mean=4200.0)},
        reason=reason)


def _reporting_task(specs=None, geom=None, rows=None, two_stage=False, plateau_tol=0.0,
                    phase2_start=None, tmpdir='.'):
    """A task carrying just enough state for the reporting hooks."""
    t = _task(dict(BASE, **(specs or {})), geom=geom)
    t.taskname = 'membrane_equilibrate'
    t.basename = os.path.join(tmpdir, 'me')
    t._mass_amu = 1.234e6
    t._two_stage = two_stage
    t._phase2_start_step = phase2_start
    t._area_plateau_tol = plateau_tol
    t._last_plateau_drift = 0.0012
    t._rows = rows if rows is not None else [(1, 5000, 5000, _joint(), 1)]
    t.register = mock.Mock()
    return t


class TestStopReasons(unittest.TestCase):

    def test_converged_prefix_is_what_the_loop_keys_off(self):
        # ChunkedEquilibrateTask.do sets self.converged from stop_reason.startswith('CONVERGED')
        r = _reporting_task()._converged_stop_reason(50_000)
        self.assertTrue(r.startswith('CONVERGED'))
        self.assertIn('density + lateral area', r)

    def test_converged_reason_reports_the_area_plateau_when_gated_on_it(self):
        r = _reporting_task(plateau_tol=0.005)._converged_stop_reason(50_000)
        self.assertIn('area plateaued', r)
        self.assertIn('quarter-drift', r)

    def test_converged_reason_omits_the_plateau_when_not_gated(self):
        self.assertNotIn('plateau', _reporting_task(plateau_tol=0.0)._converged_stop_reason(50_000))

    def test_ceiling_reason_names_both_observables(self):
        r = _reporting_task()._ceiling_stop_reason(90_000, 90_000)
        self.assertTrue(r.startswith('CEILING'))
        self.assertIn('density drift', r)
        self.assertIn('area drift', r)
        self.assertIn('may not have settled', r)

    def test_ceiling_reason_says_which_gate_was_unmet(self):
        rows = [(1, 5000, 5000,
                 _joint(converged=False, area=_report(mean=4200.0, precision_met=False)), 2)]
        r = _reporting_task(rows=rows)._ceiling_stop_reason(90_000, 90_000)
        self.assertIn('UNMET', r, 'a ceilinged run must say which precision gate failed')

    def test_ceiling_reason_survives_having_no_rows(self):
        # a run that ceilings before any chunk completes must still produce a reason
        r = _reporting_task(rows=[])._ceiling_stop_reason(0, 90_000)
        self.assertTrue(r.startswith('CEILING'))
        self.assertIn('n/a', r)

    def test_blowup_reason_and_error_mention_both_observables(self):
        t = _reporting_task()
        self.assertIn('non-finite', t._blowup_stop_reason(1234))
        self.assertIn('area', t._blowup_stop_reason(1234))
        self.assertIn('non-finite', t._blowup_error(1234))


class TestConvergenceReportFile(unittest.TestCase):
    """The ``-membrane.dat`` report is a deliverable: it is what a user inspects to decide whether
    the bilayer actually relaxed, so its header must not misdescribe how APL was computed."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = self._tmp.name
        self.addCleanup(self._tmp.cleanup)

    def _write(self, **kw):
        t = _reporting_task(tmpdir=self.tmpdir, **kw)
        t._write_report('CONVERGED at step 50000: test')
        with open(f'{t.basename}-membrane.dat') as f:
            return f.read(), t

    def test_header_flags_that_apl_is_not_protein_corrected_without_geometry(self):
        text, _t = self._write(specs=dict(lipids_per_leaflet=70))
        self.assertIn('NOT protein-corrected', text,
                      'without measured geometry the report must not imply a corrected APL')

    def test_header_describes_protein_correction_when_geometry_was_measured(self):
        geom = LeafletGeometry(n_lower=70, n_upper=70, a_prot_lower=140.0, a_prot_upper=700.0,
                               midplane_z=0.0, n_lipids_total=140)
        text, _t = self._write(geom=geom)
        self.assertIn('protein-corrected per leaflet', text)
        self.assertIn('140', text)          # the measured footprint appears
        self.assertNotIn('NOT protein-corrected', text)

    def test_two_stage_protocol_is_described_with_its_handoff_step(self):
        text, _t = self._write(two_stage=True, phase2_start=20_000)
        self.assertIn('two-stage', text)
        self.assertIn('20000', text)

    def test_two_stage_still_in_stage_one_says_so(self):
        text, _t = self._write(two_stage=True, phase2_start=None)
        self.assertIn('still in stage 1', text)

    def test_single_stage_protocol_is_described(self):
        text, _t = self._write(two_stage=False)
        self.assertIn('single-stage', text)

    def test_the_stop_reason_is_recorded_in_the_file(self):
        text, _t = self._write()
        self.assertIn('# stop: CONVERGED at step 50000: test', text)

    def test_one_data_row_per_chunk_with_the_stage_column(self):
        rows = [(1, 5000, 5000, _joint(), 1), (2, 12000, 7000, _joint(), 2)]
        text, _t = self._write(rows=rows)
        data = [l for l in text.splitlines() if not l.startswith('#')]
        self.assertEqual(len(data), 2)
        self.assertTrue(data[0].split()[0] == '1' and data[1].split()[0] == '2',
                        f'stage column wrong: {data}')

    def test_the_report_is_registered_as_an_artifact(self):
        _text, t = self._write()
        t.register.assert_called_once()
        self.assertEqual(t.register.call_args.kwargs.get('key'), 'membrane_report')


# --- _setup: protocol resolution and graceful degradation ------------------------------------
#
# _setup decides *what the task will do*: one stage or two, which observables gate convergence,
# and whether APL can be protein-corrected.  Those decisions are made from the specs and from
# whether the leaflet geometry could be measured -- and the geometry step is allowed to fail, in
# which case the report must fall back to a naive APL rather than claim a correction it did not
# make.  Everything here is resolution logic; the file reads it depends on are stubbed.


def _setup_task(specs=None, geom=None, geom_raises=None, has_coords=True):
    """Run the real ``_setup`` with its file-reading dependencies stubbed."""
    t = MembraneEquilibrateTask.__new__(MembraneEquilibrateTask)
    t.specs = dict(BASE, **(specs or {}))
    t.taskname = 'membrane_equilibrate'

    state = mock.Mock()
    state.psf.path = '/stub.psf'
    if has_coords:
        state.pdb.path = '/stub.pdb'
    else:
        state.pdb = None
        state.coor = None

    geom_mock = (mock.Mock(side_effect=geom_raises) if geom_raises
                 else mock.Mock(return_value=geom))
    with mock.patch.object(me_mod, 'total_mass_amu', return_value=1.0e6), \
         mock.patch.object(me_mod, 'total_atoms', return_value=50_000), \
         mock.patch.object(me_mod, 'membrane_leaflet_geometry', geom_mock):
        t._setup(state, min_steps=4000)
    return t


_GEOM = LeafletGeometry(n_lower=70, n_upper=70, a_prot_lower=140.0, a_prot_upper=700.0,
                        midplane_z=0.0, n_lipids_total=140)


class TestSetupProtocolResolution(unittest.TestCase):

    def test_two_stage_is_the_default_and_gates_stage_one_on_density_alone(self):
        t = _setup_task(geom=_GEOM)
        self.assertTrue(t._two_stage)
        self.assertEqual(t._phase, 1)
        self.assertEqual(set(t._conv.monitors), {'density'},
                         'stage 1 pins the area, so area is not yet a meaningful observable')

    def test_single_stage_gates_on_both_observables_at_once(self):
        t = _setup_task(dict(two_stage=False), geom=_GEOM)
        self.assertFalse(t._two_stage)
        self.assertEqual(t._phase, 2)
        self.assertEqual(set(t._conv.monitors), {'density', 'area'})

    def test_area_gets_its_own_larger_step_floor(self):
        # area is a slow, soft mode; a bigger floor stops it converging on a momentary flat spot
        t = _setup_task(dict(area_min_steps=20_000), geom=_GEOM)
        self.assertEqual(t._min_steps, 4000)
        self.assertEqual(t._area_min, 20_000)

    def test_area_floor_never_falls_below_the_task_minimum(self):
        t = _setup_task(dict(area_min_steps=100), geom=_GEOM)
        self.assertEqual(t._area_min, 4000)

    def test_constant_ratio_defaults_on_and_can_be_disabled(self):
        self.assertTrue(_setup_task(geom=_GEOM)._constant_ratio)
        self.assertFalse(_setup_task(dict(constant_ratio=False), geom=_GEOM)._constant_ratio)

    def test_plateau_gate_is_off_unless_asked_for(self):
        self.assertEqual(_setup_task(geom=_GEOM)._area_plateau_tol, 0.0)
        self.assertAlmostEqual(_setup_task(dict(area_plateau_tol=0.005),
                                           geom=_GEOM)._area_plateau_tol, 0.005)


class TestSetupGeometryDegradation(unittest.TestCase):
    """Measuring leaflet geometry is allowed to fail; claiming a protein correction that was never
    computed is not."""

    def test_measured_geometry_is_used_when_available(self):
        t = _setup_task(geom=_GEOM)
        self.assertIs(t._geom, _GEOM)

    def test_a_failed_measurement_degrades_to_no_geometry(self):
        t = _setup_task(dict(lipids_per_leaflet=70),
                        geom_raises=RuntimeError('no lipids found'))
        self.assertIsNone(t._geom, 'a failed measurement must not leave a partial geometry')
        self.assertEqual(t._lipids_per_leaflet, 70, 'the spec fallback should still be available')

    def test_missing_coordinates_degrade_rather_than_raise(self):
        t = _setup_task(dict(lipids_per_leaflet=70), has_coords=False)
        self.assertIsNone(t._geom)

    def test_degraded_geometry_reports_a_naive_apl(self):
        # the end-to-end consequence: the report header must say so, and APL must be uncorrected
        t = _setup_task(dict(lipids_per_leaflet=70),
                        geom_raises=RuntimeError('boom'))
        lo, up = t._apl_pair(4200.0)
        self.assertEqual((lo, up), (60.0, 60.0), 'without geometry APL is area/lipids_per_leaflet')
