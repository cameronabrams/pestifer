# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The calibration audit computes its own number instead of quoting the verdict it audits.

``_relaxed_area_drift`` returned the convergence gate's stored ``_last_plateau_drift`` whenever a
plateau-gated stage reported CONVERGED. The gate declares convergence **only** when that number
is below ``area_plateau_tol`` (0.005), and the downstream check asks whether it EXCEEDS
``_AREA_CONVERGENCE_TOL`` (also 0.005). So for any converged patch the check was
``|x| > 0.005`` where ``|x| < 0.005`` by construction -- a tautology.

The fatal calibration gate therefore fired only on ``CEILING``, which was already loud, and was
**structurally blind to premature convergence** -- a patch that stops early and certifies a
confident, specific, wrong APL. That is the dangerous case. Measured: with example 17's own
settings a patch converged at 88,210 steps certifying 50.53 A^2 against a truth near 43.

Recomputing cannot repair the statistic -- it is the same formula -- but it removes the property
that made the check unfalsifiable, and it is a prerequisite for ever gating on something better.

**Why the statistic itself cannot answer the question**, recorded here because it is the reason
the gate is not simply re-tuned: ``_area_plateau_drift`` is a difference of adjacent-window
MEANS, so its value scales with the length of the window. The same steady descent reads flatter
the earlier you stop -- the metric rewards stopping early. Across 28 stage-2 series measured
2026-10-14, runs the gate called CONVERGED were changing area FASTER (median 1.33 % per 100k
steps) than runs that merely hit their ceiling (0.76). The verdict is anti-correlated with the
property it names.

:meth:`MembraneEquilibrateTask.area_rate` is the independent measure. It is reported in the
convergence report and the run record and is deliberately **not** gated on: switching the gate to
a rate would mark most membrane calibrations unconverged, which is probably correct and is a
decision about what the patch protocol promises rather than a tolerance to tune.
"""
import inspect
import unittest
from unittest import mock

from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask, _as_float


class TestTheAuditDoesNotQuoteTheGate(unittest.TestCase):

    @staticmethod
    def _code(text):
        """Executable lines only.

        Three negative controls in this repository have now matched a COMMENT that quotes the
        very string they forbid, rather than the code. A guard aimed one level off its target
        passes for the wrong reason, which is the failure it exists to prevent.
        """
        return '\n'.join(ln.split('#', 1)[0] for ln in text.splitlines()
                         if ln.split('#', 1)[0].strip())

    def setUp(self):
        self.src = self._code(inspect.getsource(MakeMembraneSystemTask._relaxed_area_drift))

    def test_it_does_not_read_the_stored_verdict(self):
        """THE regression. `_last_plateau_drift` is the gate's own sub-tolerance number."""
        self.assertNotIn('_last_plateau_drift', self.src,
                         'the audit is quoting the verdict it audits again')

    def test_it_recomputes_from_the_series(self):
        self.assertIn('_area_plateau_drift', self.src)
        self.assertIn('callable', self.src,
                      'the recompute must call the method, not read an attribute of the same name')

    def test_a_converged_stage_still_yields_a_number(self):
        """POSITIVE CONTROL: breaking the shortcut must not make the audit return None and
        silently disable the (still useful) CEILING branch of the gate."""
        t = mock.Mock(_area_plateau_tol=0.005, converged=True)
        t._area_plateau_drift = lambda: -0.0042
        task = MakeMembraneSystemTask.__new__(MakeMembraneSystemTask)
        self.assertAlmostEqual(
            MakeMembraneSystemTask._relaxed_area_drift(task, [t], t, 'patchA'), -0.0042)

    def test_a_stage_that_stored_a_DIFFERENT_value_is_ignored(self):
        """The recompute must win over the stored number, which is the whole point."""
        t = mock.Mock(_area_plateau_tol=0.005, converged=True, _last_plateau_drift=0.0001)
        t._area_plateau_drift = lambda: -0.0400
        task = MakeMembraneSystemTask.__new__(MakeMembraneSystemTask)
        self.assertAlmostEqual(
            MakeMembraneSystemTask._relaxed_area_drift(task, [t], t, 'patchA'), -0.0400)


class TestTheIndependentRateIsReportedNotGated(unittest.TestCase):

    def test_the_rate_exists_and_is_normalized_per_100k_steps(self):
        from pestifer.tasks.membrane_equilibrate import MembraneEquilibrateTask
        src = inspect.getsource(MembraneEquilibrateTask.area_rate)
        self.assertIn('polyfit', src, 'a rate must come from a slope, not a window difference')
        self.assertIn('1.0e5', src)

    def test_the_rate_reaches_the_run_record(self):
        src = inspect.getsource(MakeMembraneSystemTask.build_grid_membrane_asymmetric)
        self.assertIn('apl_upper_area_rate', src)
        self.assertIn('apl_lower_area_rate', src)

    def test_the_gate_still_tests_the_DRIFT_not_the_rate(self):
        """Deliberate. Gating on the rate would fail most membrane builds -- a decision about the
        protocol, not a bug fix. Pinned so it is a choice rather than a drift."""
        src = inspect.getsource(MakeMembraneSystemTask.build_grid_membrane_asymmetric)
        gate = src[src.index('unconverged = []'):src.index('if self.embedding:')]
        self.assertIn('area_drift', gate)
        self.assertNotIn('area_rate', gate)


class TestTheRateIsAdvisoryAndCannotBreakARun(unittest.TestCase):
    """`_measure_thickness` in the same module is wrapped so an advisory observable cannot end a
    2.5M-step equilibration -- after an AttributeError from a half-initialized local escaped into
    the convergence loop. `area_rate` is the same kind of thing and needs the same property: it
    reads series that `_setup` creates, and the report writer calls it unconditionally."""

    def _bare(self):
        from pestifer.tasks.membrane_equilibrate import MembraneEquilibrateTask
        return MembraneEquilibrateTask.__new__(MembraneEquilibrateTask)

    def test_a_task_that_never_ran_a_stage_returns_None(self):
        self.assertIsNone(self._bare().area_rate())

    def test_too_few_samples_returns_None(self):
        t = self._bare()
        t._all_t, t._all_a, t._phase2_start_step = [1.0, 2.0], [10.0, 10.0], None
        self.assertIsNone(t.area_rate())

    def test_mismatched_series_lengths_return_None(self):
        """NEGATIVE CONTROL: ragged series must not raise or silently zip short."""
        t = self._bare()
        t._all_t, t._all_a, t._phase2_start_step = list(range(20)), [1.0] * 5, None
        self.assertIsNone(t.area_rate())

    def test_a_steady_decline_reports_the_right_rate(self):
        """POSITIVE CONTROL: it must measure the right number, not merely never raise.

        A gentle decline on purpose. The rate normalizes by the LATE-HALF mean area, so a series
        that loses a large fraction of its area makes the expected value depend on which window
        normalizes it -- and a first draft of this test asserted -100 %/100k for a series whose
        true answer was -141.8, because it divided by the starting area instead. Kept gentle so
        the expectation is derivable without reproducing the implementation.
        """
        t = self._bare()
        n = 40
        slope = -1.0e-5                       # area units per step
        t._all_t = [float(i * 1000) for i in range(n)]
        t._all_a = [100.0 + slope * (i * 1000) for i in range(n)]      # 100.0 -> 99.61
        t._phase2_start_step = None
        r = t.area_rate()
        self.assertIsNotNone(r)
        self.assertLess(r, 0.0)
        # 1e5 * slope / ~100 area units, as a percent -> about -1.0 %/100k
        self.assertAlmostEqual(r, -1.0, delta=0.05)

    def test_a_flat_series_reports_essentially_zero(self):
        """The other end: a settled area must not read as drifting."""
        t = self._bare()
        t._all_t = [float(i * 1000) for i in range(40)]
        t._all_a = [100.0] * 40
        t._phase2_start_step = None
        self.assertAlmostEqual(t.area_rate(), 0.0, places=6)


class TestAsFloatRejectsStandIns(unittest.TestCase):
    """A stand-in answers `getattr` with something that is neither None nor a number, and
    formatting or `float()`-ing it raises from inside a build. `_divisor` in this same module
    already carries that lesson; this is the second site to need it."""

    def test_a_mock_is_not_a_number(self):
        self.assertIsNone(_as_float(mock.Mock()))

    def test_a_bool_is_not_a_number_either(self):
        self.assertIsNone(_as_float(True))

    def test_real_numbers_pass(self):
        self.assertEqual(_as_float(2), 2.0)
        self.assertEqual(_as_float(-0.5), -0.5)
        self.assertIsNone(_as_float(None))


if __name__ == '__main__':
    unittest.main()
