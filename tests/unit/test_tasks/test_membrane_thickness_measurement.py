# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
P-P thickness is measured per chunk by ``MembraneEquilibrateTask`` and reported, but **does not
gate**.  These tests pin both halves of that: that it is actually measured and reaches the log and
the report, and that nothing about it can change whether a build certifies or fail a build when it
cannot be measured.

The second half is the one that matters while the observable is advisory.  A measurement bolted
onto the convergence loop that can raise is a measurement that can take down a 2.5M-step
equilibration for no scientific reason.
"""
import contextlib
import os
import tempfile
import unittest

from unittest.mock import MagicMock, patch

from pestifer.tasks.membrane_equilibrate import MembraneEquilibrateTask
from pestifer.util.density_convergence import MembraneThickness


@contextlib.contextmanager
def _in_tmpdir():
    """A scratch cwd that is always restored.

    conftest.py gives each test module its own working directory; a test that chdir'd away and
    did not come back would silently relocate every test that ran after it in the same process.
    """
    prev = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        try:
            os.chdir(d)
            yield d
        finally:
            os.chdir(prev)


def _task(**attrs):
    """A task with only the attributes `_measure_thickness` / `_thickness_note` touch."""
    t = MembraneEquilibrateTask.__new__(MembraneEquilibrateTask)
    t.taskname = 'membrane_equilibrate'
    t.basename = 'me'
    t._psf_path = 'system.psf'
    t._all_pp = []
    t._last_pp = None
    t._mon_pp = MagicMock()
    for k, v in attrs.items():
        setattr(t, k, v)
    return t


def _result(thickness=38.5, note='', n_lower=64, n_upper=64, n_unresolved=0):
    return MembraneThickness(thickness, None, None, n_upper, n_lower, 0, n_unresolved, 0.0, note)


class TestThicknessIsMeasured(unittest.TestCase):

    def test_a_measured_thickness_is_recorded_and_monitored(self):
        t = _task()
        with _in_tmpdir():
            open('me.coor', 'wb').close()
            with patch('pestifer.tasks.membrane_equilibrate.membrane_pp_thickness',
                       return_value=_result(38.5)):
                t._measure_thickness(120000)
        self.assertEqual(t._all_pp, [(120000, 38.5)])
        t._mon_pp.add_samples.assert_called_once_with([120000.0], [38.5])

    def test_it_reaches_the_log_line(self):
        t = _task(_all_pp=[(1, 38.5)] * 2, _last_pp=_result(38.5))
        note = t._thickness_note()
        self.assertIn('P-P=38.50 A', note)
        self.assertIn('64/64 P', note)

    def test_ambiguous_lipids_are_named_in_the_log_not_hidden(self):
        """Excluding a cardiolipin silently would make the number look better than it is."""
        t = _task(_all_pp=[(1, 38.5)], _last_pp=_result(38.5, n_unresolved=3))
        self.assertIn('3 lipid(s) with ambiguous phosphorus excluded', t._thickness_note())

    def test_an_unmeasurable_thickness_says_nothing_rather_than_zero(self):
        t = _task(_last_pp=_result(None, note='no resolved phosphorus in upper'))
        self.assertEqual(t._thickness_note(), '')


class TestThicknessCannotBreakABuild(unittest.TestCase):
    """NEGATIVE CONTROLS.  While this observable is advisory, every failure of it must be inert."""

    def test_a_raising_measurement_is_swallowed(self):
        t = _task()
        with _in_tmpdir():
            open('me.coor', 'wb').close()
            with patch('pestifer.tasks.membrane_equilibrate.membrane_pp_thickness',
                       side_effect=ValueError('not a membrane system')):
                t._measure_thickness(1000)       # must not raise
        self.assertEqual(t._all_pp, [])
        self.assertIsNone(t._last_pp)

    def test_a_missing_coor_is_not_fatal(self):
        t = _task()
        with _in_tmpdir():
            t._measure_thickness(1000)           # no me.coor written at all
        self.assertEqual(t._all_pp, [])

    def test_an_undefined_thickness_is_not_recorded_as_a_number(self):
        t = _task()
        with _in_tmpdir():
            open('me.coor', 'wb').close()
            with patch('pestifer.tasks.membrane_equilibrate.membrane_pp_thickness',
                       return_value=_result(None, note='no resolved phosphorus in both leaflets')):
                t._measure_thickness(1000)
        self.assertEqual(t._all_pp, [], 'an undefined thickness must not enter the series')
        t._mon_pp.add_samples.assert_not_called()

    def test_thickness_is_absent_from_the_convergence_gate(self):
        """The whole point of landing this measure-only.

        If thickness ever appears in JointConvergence, certification behaviour -- and the runtime
        of every membrane build, including ex16's, whose gate is already correct -- changes.  That
        is a deliberate second step, so this test fails if it happens by accident.
        """
        import inspect
        src = inspect.getsource(MembraneEquilibrateTask)
        joint = [ln for ln in src.splitlines() if 'JointConvergence(' in ln]
        self.assertTrue(joint, 'POSITIVE CONTROL: no JointConvergence construction found')
        for ln in joint:
            self.assertNotIn('_mon_pp', ln)
            self.assertNotIn('thickness', ln.lower())

    def test_the_monitor_exists_so_gating_is_a_one_line_change(self):
        """Counterpart to the test above: measure-only must not mean 'no statistics'.

        The thickness series runs through a real DensityConvergenceMonitor -- autocorrelation
        corrected, upper-confidence-bounded -- so switching it on later is adding it to
        JointConvergence, not writing a new criterion under time pressure.
        """
        import inspect
        src = inspect.getsource(MembraneEquilibrateTask._setup)
        self.assertIn('self._mon_pp = DensityConvergenceMonitor', src)


if __name__ == '__main__':
    unittest.main()
