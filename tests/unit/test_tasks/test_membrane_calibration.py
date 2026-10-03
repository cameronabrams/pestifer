# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Two bugs reported 2026-10-03 against 3.24.1, from an asymmetric pseudo-C3 Env build on picotte
(job 26186243).  Both live in ``make_membrane_system`` and both are silent in their own way.

**The diagnostic re-solvated an already-solvated quilt.**  ``equilibrate_bilayer`` appended a
``solvate`` task unconditionally -- right for a freshly gridded patch, whose chambers the packer
leaves empty, and wrong for the finished quilt the differential-stress diagnostic hands back to
it.  psfgen stopped with ``duplicate segment key WT1`` after 19 h of CPU, so the diagnostic could
never produce the leaflet-stress number it exists for.

**The calibrated APL divided by the requested lipid count.**  ``ring_check`` deletes piercing
lipids during the patch protocol, so a relaxed patch holds fewer than ``patch_nlipids``.  Dividing
its area by the request gives an APL up to ~7% low, and the stress-free counts derived from it
over-fill that leaflet: 744/691 lipids placed where the measured counts give 737/660, leaving the
leaflets mismatched by ~3.7%.  That is a built-in differential stress -- the quantity the
calibration exists to remove -- and nothing in the build reports it.
"""
import unittest

from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask


def _task(**attrs):
    t = MakeMembraneSystemTask.__new__(MakeMembraneSystemTask)
    t.taskname = 'make_membrane_system'
    t.basename = 'mms'
    t.embedding = False
    t.bilayer_specs = dict(patch_nlipids=dict(upper=100, lower=100), npatch=[1, 1])
    for k, v in attrs.items():
        setattr(t, k, v)
    return t


def _patch(area, n_leaflet=None, area_drift=0.0):
    p = SimpleNamespace(area=area, area_drift=area_drift)
    if n_leaflet is not None:
        p.n_leaflet = n_leaflet
    return p


class TestSolvateIsSkippedWhenAlreadySolvated(unittest.TestCase):
    """Guarded on the STRUCTURE, not on which caller it is."""

    def _tasklist(self, has_water):
        """The task list `equilibrate_bilayer` would build, without running anything."""
        import inspect
        src = inspect.getsource(MakeMembraneSystemTask.equilibrate_bilayer)
        self.assertIn('psf_contains_water', src,
                      'the solvate step is no longer gated on whether the state holds water')
        # the gate must wrap the solvate append, not merely be mentioned
        gate = src.index('psf_contains_water')
        solv = src.index("{'solvate'")
        self.assertLess(gate, solv, 'psf_contains_water is checked after solvate is appended')
        return src

    def test_the_solvate_step_is_gated(self):
        self._tasklist(has_water=True)

    def test_the_gate_reads_the_psf_not_the_caller(self):
        """A flag passed down by the diagnostic would fix that one call site and miss the next.

        NEGATIVE CONTROL for the shape of the original bug: the decision has to come from the
        structure, so any future path handing over a solvated state is covered without being
        enumerated.
        """
        import inspect
        src = inspect.getsource(MakeMembraneSystemTask.equilibrate_bilayer)
        line = next(ln for ln in src.splitlines() if 'psf_contains_water' in ln and 'if' in ln)
        self.assertIn('state.psf', line,
                      f'the gate does not read the incoming PSF: {line.strip()!r}')


class TestCalibratedAPLUsesTheMeasuredCount(unittest.TestCase):

    def _apls(self, t):
        """Run just the divisor/APL arithmetic out of build_grid_membrane_asymmetric."""
        captured = {}

        def fake_grid(counts, box_SAPL, aspect):
            captured['counts'] = counts

        t._npatch = lambda: [1, 1]
        t._grid_membrane = fake_grid
        t.specs = {}
        with patch.object(MakeMembraneSystemTask, 'diagnose_differential_stress'):
            MakeMembraneSystemTask.build_grid_membrane_asymmetric(t)
        return captured['counts']

    def test_each_patch_divides_by_its_own_measured_count(self):
        """patchA ended with 99 lipids/leaflet and patchB with 95.5, not the 100 requested."""
        t = _task(patchA=_patch(area=4134.0, n_leaflet=99.0),
                  patchB=_patch(area=4451.0, n_leaflet=95.5))
        counts = self._apls(t)
        # APL_upper = 4134/99 = 41.76, APL_lower = 4451/95.5 = 46.61 (the reported values)
        # n_lower/n_upper must follow apl_upper/apl_lower
        ratio = counts['lower'] / counts['upper']
        self.assertAlmostEqual(ratio, (4134.0 / 99.0) / (4451.0 / 95.5), places=2)

    def test_dividing_by_the_request_overfills_the_leaflet(self):
        """POSITIVE CONTROL: the old arithmetic really does bias the counts.

        Without it, the test above could pass against a formula that ignored the counts entirely.
        """
        measured = self._apls(_task(patchA=_patch(4134.0, n_leaflet=99.0),
                                    patchB=_patch(4451.0, n_leaflet=95.5)))
        requested = self._apls(_task(patchA=_patch(4134.0, n_leaflet=100.0),
                                     patchB=_patch(4451.0, n_leaflet=100.0)))
        self.assertGreater(requested['lower'], measured['lower'],
                           'dividing by the requested count must over-fill the lower leaflet')

    def test_the_lower_patch_is_not_divided_by_the_upper_request(self):
        """The second, unreported half of the bug.

        `n_cal = patch_nlipids['upper']` was the divisor for BOTH patches, so a config asking for
        different upper and lower counts mis-scaled the lower leaflet too.  The reporter's run had
        both at 100 and could not have seen it.
        """
        t = _task(patchA=_patch(4134.0), patchB=_patch(4451.0))   # no measured counts -> fallback
        t.bilayer_specs = dict(patch_nlipids=dict(upper=100, lower=80), npatch=[1, 1])
        counts = self._apls(t)
        expected = (4134.0 / 100.0) / (4451.0 / 80.0)     # each by ITS OWN request
        wrong = (4134.0 / 100.0) / (4451.0 / 100.0)       # both by the upper request
        self.assertAlmostEqual(counts['lower'] / counts['upper'], expected, places=2)
        self.assertNotAlmostEqual(counts['lower'] / counts['upper'], wrong, places=2)

    def test_box_sizing_still_uses_the_requested_count(self):
        """What must NOT change: `patch_nlipids x npatch` is the system the user asked for.

        Only the APL divisor cares what ring_check deleted; sizing the box from the measured count
        would quietly hand back a smaller membrane than was requested.
        """
        import inspect
        src = inspect.getsource(MakeMembraneSystemTask.build_grid_membrane_asymmetric)
        line = next(ln for ln in src.splitlines()
                    if 'npatch[0] * npatch[1]' in ln and 'n_upper' in ln)
        self.assertIn("patch_nlipids['upper']", line,
                      f'box sizing no longer uses the requested count: {line.strip()!r}')


if __name__ == '__main__':
    unittest.main()
