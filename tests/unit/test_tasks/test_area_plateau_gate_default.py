# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The area-plateau gate is ON by default for a calibration patch, and OFF for the quilt.

A calibration patch's APL sizes the whole grid. Left ungated, ``membrane_equilibrate`` could
certify that stage on a small *local* slope at ``area_drift_tol`` (0.012) -- looser than the
0.005 ``make_membrane_system`` later applies to the same system when judging whether the
calibrated APL is trustworthy. Two tolerances for one judgment, and only the looser one enforced.

Measured consequence, reported from salsa's fleet 2026-10-08: a patch certified an APL of 45.62
and then condensed to 43.47 -- **4.7% early**, roughly nine times the downstream tolerance -- and
the leaflet sized from it came out ~5% under-filled (674 lipids placed against ~708).

``area_plateau_tol`` closes it by requiring the *cumulative* quarter-over-quarter drift of the
stage-2 area series below the same 0.005, which is the exact metric the downstream check uses, so
the two now agree by construction rather than by luck.

Why ON for the patch and OFF for the quilt: nothing downstream is sized from the quilt's area, so
a slightly early stop there costs a little extra relaxation later rather than a mis-sized system.
The gate is strictly stricter -- it can only ever *refuse* a convergence the local gate would have
granted -- so defaulting it on cannot introduce a new wrong answer, only more steps.

``max_steps`` had to move with it. It is a CEILING on a loop that stops itself at convergence, so
a fast patch pays nothing for headroom. ex17's calibration patch across three replicas took
236870 and 372360 steps and was still descending at 1200000, so at the old 500000 the new gate
would have produced a CEILING warning instead of a convergence on the exact systems it exists
for.

No shipped example changes behavior: ex16 and ex17 are the only configs with relaxation
protocols and both already set ``area_plateau_tol`` explicitly.
"""
import unittest
from importlib.resources import files as pkg_files

import yaml

_SPEC = yaml.safe_load(open(str(pkg_files('pestifer.schema').joinpath('base.yaml'))))

#: the downstream reliability tolerance the gate has to agree with
_DOWNSTREAM_TOL = 0.005


def _find(node, name, _seen=None):
    """Every declaration of `name`, in document order."""
    out = []
    if isinstance(node, dict):
        if node.get('name') == name:
            out.append(node)
        for v in node.values():
            out += _find(v, name)
    elif isinstance(node, list):
        for v in node:
            out += _find(v, name)
    return out


def _stage(name):
    """The relaxation_protocols stage declaration called `name`."""
    for cand in _find(_SPEC, name):
        if isinstance(cand.get('attributes'), list) or isinstance(cand.get('default'), list):
            return cand
    raise AssertionError(f'no {name} stage in the schema')


def _attr_of_stage(stage_name, attr):
    hits = _find(_stage(stage_name), attr)
    assert hits, f'{stage_name} has no {attr}'
    return hits[0]


class TestTheDefaults(unittest.TestCase):

    def test_the_patch_stage_gates_on_the_plateau_by_default(self):
        self.assertEqual(_attr_of_stage('patch', 'area_plateau_tol')['default'], _DOWNSTREAM_TOL)

    def test_the_quilt_stage_does_not(self):
        """Deliberate, not an oversight -- nothing downstream is sized from the quilt's area."""
        self.assertEqual(_attr_of_stage('quilt', 'area_plateau_tol')['default'], 0.0)

    def test_the_gate_matches_the_tolerance_that_later_judges_the_result(self):
        """The whole point. If these two drift apart, the stage certifies at one tolerance and is
        graded at another -- which is the defect this closes."""
        from pestifer.tasks.make_membrane_system import _AREA_CONVERGENCE_TOL
        self.assertEqual(_attr_of_stage('patch', 'area_plateau_tol')['default'],
                         _AREA_CONVERGENCE_TOL)

    def test_the_patch_ceiling_was_raised_to_match(self):
        """At 500000 the new gate would report CEILING on the cholesterol-rich patches it is for."""
        self.assertGreaterEqual(_attr_of_stage('patch', 'max_steps')['default'], 1_200_000)

    def test_the_quilt_ceiling_was_NOT_raised(self):
        """NEGATIVE CONTROL for a blanket edit: the quilt's gate stays off, so its budget stays."""
        self.assertEqual(_attr_of_stage('quilt', 'max_steps')['default'], 500000)

    def test_the_local_tolerance_is_still_looser_than_the_gate(self):
        """Recorded because it is the asymmetry that caused this, and it is NOT removed -- the
        local gate still certifies at 0.012; it simply can no longer certify alone."""
        self.assertGreater(_attr_of_stage('patch', 'area_drift_tol')['default'], _DOWNSTREAM_TOL)

    def test_every_declaration_carries_text(self):
        """A schema key with no text: crashes the RTD config-ref build."""
        for d in _find(_SPEC, 'area_plateau_tol'):
            self.assertIn('text', d)
            self.assertEqual(d['type'], 'float')


class TestTheGateIsStrictlyStricter(unittest.TestCase):
    """It can only ever REFUSE a convergence the local gates granted.

    That is what makes defaulting it on safe: it cannot produce a new class of wrong answer, only
    more steps. If it ever moves outside the `if converged` branch, it could start *granting*
    convergence, and this stops being true.
    """

    def test_the_gate_only_runs_when_the_local_gates_already_passed(self):
        import inspect
        from pestifer.tasks.membrane_equilibrate import MembraneEquilibrateTask
        src = inspect.getsource(MembraneEquilibrateTask)
        i = src.index('_area_plateau_tol > 0.0')
        window = src[max(0, i - 400):i]
        self.assertIn('converged = jr.converged', window)
        self.assertIn('if converged and', window,
                      'the plateau gate is no longer guarded by the local gates having passed, '
                      'so it could now grant convergence rather than only withhold it')


if __name__ == '__main__':
    unittest.main()
