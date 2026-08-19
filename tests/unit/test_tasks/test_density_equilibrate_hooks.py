# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Unit tests for the density-specific hooks of ``density_equilibrate``.

``ChunkedEquilibrateTask`` supplies the loop (tested in ``test_equilibrate_loop.py``) and
``DensityConvergenceMonitor`` supplies the criterion (tested in ``test_density_convergence.py``).
What sits between them -- reading a chunk's ``.xst``, turning cell volumes into densities, feeding
the monitor, and turning its verdict into the stop reason a user reads and ``run-record.json``
stores -- had no test, because it only ran under a live NAMD.

It does not need one.  A ``.xst`` is a text file, so a density trajectory can simply be written
and handed to the hook.  These tests therefore drive real convergence decisions from synthetic
series: a settled box, a still-densifying box, and a box that blew up.
"""

import math
import os
import tempfile
import unittest
from unittest import mock

import numpy as np

from pestifer.tasks import density_equilibrate as DE
from pestifer.tasks.density_equilibrate import DensityEquilibrateTask


_XST_HEADER = ('# NAMD extended system trajectory file\n'
               '#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z\n')

def _write_xst(path, steps, volumes):
    """Write a NAMD .xst whose cells have the given volumes (cubic cells, so V = L**3)."""
    with open(path, 'w') as f:
        f.write(_XST_HEADER)
        for s, v in zip(steps, volumes):
            L = float(v) ** (1.0 / 3.0) if math.isfinite(v) and v > 0 else float('nan')
            f.write(f'{s:d} {L:.6f} 0 0 0 {L:.6f} 0 0 0 {L:.6f} 0 0 0\n')


def _task(tmpdir, **spec_overrides):
    """A DensityEquilibrateTask with its hooks live but its file/pipeline deps stubbed."""
    t = DensityEquilibrateTask.__new__(DensityEquilibrateTask)
    t.taskname = 'density_equilibrate'
    t.basename = os.path.join(tmpdir, 'de')
    t.specs = dict(drift_tol=2e-3, precision_p=3.0, drift_conf=0.0, burn_in=0,
                   window_frac=1.0, n_consecutive=3, autocorr_reliability=6.0,
                   **spec_overrides)
    t._rows = []
    state = mock.Mock()
    state.psf.path = os.path.join(tmpdir, 'stub.psf')
    with mock.patch.object(DE, 'total_mass_amu', return_value=1.0e6), \
         mock.patch.object(DE, 'total_atoms', return_value=10_000):
        t._setup(state, min_steps=0)
    return t


def _feed(task, tmpdir, volumes, step0=0, stride=100, chunk=1):
    """Write one chunk's .xst and push it through ``_ingest_chunk``; return (blowup, converged)."""
    xst = os.path.join(tmpdir, f'chunk{chunk}.xst')
    steps = [step0 + i * stride for i in range(len(volumes))]
    _write_xst(xst, steps, volumes)
    return task._ingest_chunk(xst, steps[-1] if steps else step0, len(volumes) * stride, chunk)


class _HookTestCase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmpdir = self._tmp.name
        self.addCleanup(self._tmp.cleanup)


class TestDensityConversion(_HookTestCase):

    def test_density_is_computed_from_the_cell_volume_and_system_mass(self):
        t = _task(self.tmpdir)
        # one flat chunk; the monitor's mean density should match the analytic value
        v = 125_000.0                      # a 50 A cube
        _feed(t, self.tmpdir, [v] * 40)
        report = t._rows[-1][3]
        expected = DE.volume_to_density(np.array([v]), t._mass_amu)[0]
        self.assertAlmostEqual(report.mean_density, float(expected), places=6)

    def test_a_chunk_with_no_frames_is_tolerated(self):
        """A short chunk can legitimately produce an .xst with no data rows yet."""
        t = _task(self.tmpdir)
        blowup, converged = _feed(t, self.tmpdir, [])
        self.assertFalse(blowup)
        self.assertFalse(converged)
        self.assertEqual(len(t._rows), 1, 'the chunk should still be recorded')


class TestConvergenceVerdicts(_HookTestCase):

    def test_a_settled_box_eventually_converges(self):
        t = _task(self.tmpdir)
        rng = np.random.default_rng(0)
        converged = False
        step = 0
        for chunk in range(1, 12):
            # a plateau: small independent noise about a fixed volume
            vols = 125_000.0 * (1.0 + rng.normal(0, 2e-5, 60))
            _blow, converged = _feed(t, self.tmpdir, vols, step0=step, chunk=chunk)
            step += 60 * 100
            if converged:
                break
        self.assertTrue(converged, f'a flat series never converged; last: {t._rows[-1][3].reason}')

    def test_a_densifying_box_does_not_converge(self):
        t = _task(self.tmpdir)
        step, converged = 0, False
        for chunk in range(1, 12):
            # a steady 1%-per-chunk contraction: far above the 0.2% drift tolerance
            base = 125_000.0 * (0.99 ** chunk)
            vols = np.linspace(base, base * 0.99, 60)
            _blow, converged = _feed(t, self.tmpdir, vols, step0=step, chunk=chunk)
            step += 60 * 100
            self.assertFalse(converged, f'a drifting series converged at chunk {chunk}')

    def test_a_non_finite_density_is_a_blowup(self):
        t = _task(self.tmpdir)
        blowup, converged = _feed(t, self.tmpdir, [125_000.0] * 20 + [float('nan')])
        self.assertTrue(blowup, 'a non-finite density must be reported as a blowup')
        self.assertFalse(converged)


class TestStopReasons(_HookTestCase):
    """The stop reason is what a user reads and what run-record.json stores, so its content is
    part of the contract, not decoration."""

    def _settled(self):
        t = _task(self.tmpdir)
        rng = np.random.default_rng(1)
        step = 0
        for chunk in range(1, 12):
            vols = 125_000.0 * (1.0 + rng.normal(0, 2e-5, 60))
            _blow, conv = _feed(t, self.tmpdir, vols, step0=step, chunk=chunk)
            step += 60 * 100
            if conv:
                return t, step
        self.fail('fixture did not converge')

    def test_converged_reason_names_the_drift_and_the_pass_count(self):
        t, step = self._settled()
        reason = t._converged_stop_reason(step)
        self.assertTrue(reason.startswith('CONVERGED'),
                        'the loop keys `converged` off this prefix')
        self.assertIn('drift', reason)
        self.assertIn('consecutive', reason)

    def test_ceiling_reason_reports_whether_the_precision_gate_was_met(self):
        t = _task(self.tmpdir)
        _feed(t, self.tmpdir, np.linspace(125_000.0, 120_000.0, 60))
        reason = t._ceiling_stop_reason(6000, 6000)
        self.assertTrue(reason.startswith('CEILING'))
        self.assertIn('UNMET', reason,
                      'a run that ceilinged on a drifting series must say the gate was unmet')
        self.assertIn('may not have settled', reason)

    def test_ceiling_reason_survives_having_no_rows(self):
        # a run that ceilings before any chunk completes must still produce a reason
        t = _task(self.tmpdir)
        self.assertIn('CEILING', t._ceiling_stop_reason(0, 6000))

    def test_blowup_reason_and_error_name_the_cause(self):
        t = _task(self.tmpdir)
        self.assertIn('blowup', t._blowup_stop_reason(1234).lower())
        self.assertIn('non-finite', t._blowup_error(1234).lower())
