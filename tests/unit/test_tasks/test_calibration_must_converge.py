# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""An unconverged calibration stops the build instead of sizing the grid from a bad number.

The whole membrane grid is sized from two numbers: each calibration patch's preferred area per
lipid. A patch that never equilibrated does not yield a slightly worse membrane — it yields one
built to the wrong size, and the failure surfaces far from its cause.

Observed 2026-10-12: an ex17 calibration patch hit its 2,500,000-step ceiling at area drift
**-2.50%**, five times the 0.005 tolerance, with precision UNMET. The build logged a warning,
sized the grid from that APL anyway, and died 2,100 steps into the quilt melt as a RATTLE
constraint failure on an arbitrary PE headgroup carbon — with nothing at the crash site pointing
back at the calibration. Two of three replicas died that way; the survivor was the one whose
patch converged quickly.

This is the judgment already made for ``NAMDScripter._verify_params_cover_psf``, and for the same
stated reason: the condition is fatal to the build regardless, so the only question is whether it
fails somewhere it can be diagnosed. That check also reports *every* failing term rather than the
first, because NAMD dying on one and hiding the rest is the thing it exists to prevent. Same here
— a build with both leaflets unconverged should say so once, not twice over two runs.
"""
import inspect
import unittest
from unittest import mock

from pestifer.core.errors import PestiferBuildError
from pestifer.tasks.make_membrane_system import (_AREA_CONVERGENCE_TOL, MakeMembraneSystemTask)


class _Patch:
    def __init__(self, drift):
        self.area_drift = drift


def _check(drift_upper, drift_lower):
    """Run only the convergence gate, with the surrounding sizing logic stubbed out."""
    src = inspect.getsource(MakeMembraneSystemTask.build_grid_membrane_asymmetric)
    start = src.index('unconverged = []')
    # Stop at the override call, not at `if self.embedding:` -- the override is downstream of
    # the gate and needs `bs`, which this harness deliberately does not provide. Extracting past
    # it would make the harness fail for a reason that has nothing to do with the gate.
    end = src.rindex('\n', 0, src.index('self._apply_leaflet_apl_override')) + 1
    body = '\n'.join(line[8:] if line.startswith('        ') else line
                     for line in src[start:end].rstrip().splitlines())
    ns = {
        'self': mock.Mock(patchA=_Patch(drift_upper), patchB=_Patch(drift_lower)),
        'apl_upper': 41.58, 'apl_lower': 43.47,
        '_AREA_CONVERGENCE_TOL': _AREA_CONVERGENCE_TOL,
        'PestiferBuildError': PestiferBuildError,
        'sA2_': 'A^2',
        'logger': mock.Mock(),
    }
    exec(compile(body, '<gate>', 'exec'), ns)
    return ns


class TestTheGate(unittest.TestCase):

    def test_two_converged_patches_pass(self):
        _check(0.001, -0.002)          # must not raise

    def test_a_patch_at_the_tolerance_passes(self):
        """Strictly greater-than, so the tolerance itself is not a failure."""
        _check(_AREA_CONVERGENCE_TOL, 0.0)

    def test_an_unconverged_upper_leaflet_stops_the_build(self):
        with self.assertRaises(PestiferBuildError) as cm:
            _check(-0.025, 0.001)
        msg = str(cm.exception)
        self.assertIn('upper', msg)
        self.assertIn('-2.50%', msg, 'the drift that failed is not in the message')
        self.assertIn('0.50%', msg, 'the tolerance it failed against is not in the message')

    def test_an_unconverged_lower_leaflet_stops_the_build(self):
        with self.assertRaises(PestiferBuildError):
            _check(0.001, 0.02)

    def test_BOTH_failing_leaflets_are_reported_at_once(self):
        """The `_verify_params_cover_psf` precedent: report them all, not the first.

        A build that fails, gets one leaflet fixed, and fails again on the other has cost two
        multi-hour runs to learn something one message could have said.
        """
        with self.assertRaises(PestiferBuildError) as cm:
            _check(-0.025, 0.02)
        msg = str(cm.exception)
        self.assertIn('upper', msg)
        self.assertIn('lower', msg)

    def test_the_message_names_the_remedy_and_the_downstream_symptom(self):
        """So the next person does not debug it as a RATTLE problem."""
        with self.assertRaises(PestiferBuildError) as cm:
            _check(-0.025, 0.0)
        msg = str(cm.exception)
        self.assertIn('relaxation_protocols.patch', msg)
        self.assertIn('area_plateau_tol', msg)
        self.assertIn('RATTLE', msg)

    def test_an_UNMEASURABLE_drift_is_not_fatal(self):
        """`None` means the metric could not be taken, which is a different condition from
        'taken and bad'. Making it fatal would break builds that never had the measurement."""
        ns = _check(None, None)          # must not raise
        self.assertTrue(ns['logger'].warning.called,
                        'an unmeasurable drift must at least be reported')


class TestItIsWiredIntoTheBuild(unittest.TestCase):
    """A gate nobody reaches is the original bug in a different costume."""

    def setUp(self):
        self.src = inspect.getsource(MakeMembraneSystemTask.build_grid_membrane_asymmetric)

    def test_the_gate_runs_before_the_grid_is_sized(self):
        i_gate = self.src.index('unconverged')
        i_grid = self.src.index('self._grid_membrane(')
        self.assertLess(i_gate, i_grid,
                        'the membrane is gridded before the calibration is checked')

    def test_it_raises_rather_than_warning(self):
        """NEGATIVE CONTROL for a revert to the old behavior."""
        gate = self.src[self.src.index('unconverged = []'):self.src.index('if self.embedding:')]
        self.assertIn('raise PestiferBuildError', gate)


if __name__ == '__main__':
    unittest.main()
