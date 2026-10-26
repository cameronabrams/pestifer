# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Two defects in ``diagnose_differential_stress``, reported 2026-10-07 against 3.25.1 from a
pseudo-asymmetric PIP2 build, and both confirmed in the code before being fixed.

**The pass ran at the wrong temperature.**  It builds its own protocol out of ``md`` stages and
named no temperature, so it took the ``md`` task's schema default.  The sting is that the two
relaxation task types disagree -- ``md`` defaults to 300 K, ``membrane_equilibrate`` to 310 K --
so a build that set no temperature anywhere still equilibrated at 310 K and was diagnosed at
300 K, with nothing in the config to show for it; and a build that set 310 K on every stage, as
shipped example 17 does, was diagnosed at 300 K anyway because the protocol never read it.  Some
of a temperature mismatch cancels in a *difference* of leaflet tensions, but only insofar as the
two leaflets share a thermal expansivity -- which is precisely what a compositionally asymmetric
bilayer does not do.  So it perturbed the one observable the pass exists to report.

**The pass advanced the membrane it was measuring.**  It runs through ``equilibrate_bilayer``,
which ends by rekeying its final state over ``quilt_state`` -- and ``quilt_state`` is what the
embed step consumes.  So the structure that went forward was the last frame of a fixed-length,
*ungated* run whose cell had also relaxed (``surfaceTensionTarget = 0``), in place of the
convergence-gated membrane.  Both the method docstring and the schema said "the membrane is not
rebuilt": true of rebuilding, and silent about the state that actually reaches the next task.
Two copies of one claim, neither re-derived.
"""
import inspect
import unittest

import yaml

from pestifer.tasks.make_membrane_system import (MakeMembraneSystemTask,
                                                 _MD_DEFAULT_TEMPERATURE,
                                                 _MEMBRANE_EQUILIBRATE_DEFAULT_TEMPERATURE)

_T = MakeMembraneSystemTask._relaxation_temperature


class TestTemperatureIsInherited(unittest.TestCase):

    def test_an_explicit_stage_temperature_wins(self):
        proto = [{'md': dict(ensemble='minimize')},
                 {'membrane_equilibrate': dict(temperature=293.0)}]
        self.assertEqual(_T(proto), 293.0)

    def test_the_last_dynamics_stage_decides(self):
        """A protocol ends where the membrane ended; earlier stages are not the state point."""
        proto = [{'md': dict(ensemble='NVT', temperature=410.0)},
                 {'md': dict(ensemble='NPT', temperature=310.0)}]
        self.assertEqual(_T(proto), 310.0)

    def test_a_trailing_minimize_is_not_a_state_point(self):
        """NEGATIVE CONTROL: a minimization has no thermostat, so it must not answer the question.

        Taking the last stage unconditionally would read a `minimize` with no temperature and
        fall through to a task default -- silently discarding the 310 K the stage before it
        actually ran at.
        """
        proto = [{'md': dict(ensemble='NPT', temperature=310.0)},
                 {'md': dict(ensemble='minimize', nsteps=500)}]
        self.assertEqual(_T(proto), 310.0)

    def test_a_defaulted_membrane_equilibrate_reports_ITS_default_not_mds(self):
        """The whole bug in one assertion.

        A quilt relaxed by a `membrane_equilibrate` that names no temperature ran at 310 K.
        Reporting 300 K here -- the `md` default -- is exactly the mismatch that was shipped.
        """
        self.assertEqual(_T([{'membrane_equilibrate': dict(max_steps=800000)}]),
                         _MEMBRANE_EQUILIBRATE_DEFAULT_TEMPERATURE)
        self.assertNotEqual(_MEMBRANE_EQUILIBRATE_DEFAULT_TEMPERATURE, _MD_DEFAULT_TEMPERATURE,
                            'if these two were equal the defect could not have existed, and this '
                            'whole test file would be vacuous')

    def test_a_defaulted_md_stage_reports_the_md_default(self):
        self.assertEqual(_T([{'md': dict(ensemble='NPT', nsteps=1000)}]), _MD_DEFAULT_TEMPERATURE)

    def test_a_protocol_with_no_dynamics_is_unanswerable(self):
        self.assertIsNone(_T([{'md': dict(ensemble='minimize')}]))
        self.assertIsNone(_T([]))
        self.assertIsNone(_T(None))
        self.assertIsNone(_T([{'ring_check': dict(segtypes=['lipid'])}]))


class TestTheMirroredDefaultsMatchTheSchema(unittest.TestCase):
    """The two constants are copies of schema values, so they can drift.  This is the guard.

    A copy of a number in another file is a claim about that file, and claims rot.  If someone
    retunes `membrane_equilibrate`'s default to 303 K, the diagnostic must not keep sampling at
    310 K because a constant here still says so.
    """

    @staticmethod
    def _default_for(task_name):
        from pestifer.core.resourcemanager import ResourceManager
        path = ResourceManager().schema_path if hasattr(ResourceManager(), 'schema_path') else None
        import pestifer, os
        path = path or os.path.join(os.path.dirname(pestifer.__file__), 'schema', 'base.yaml')
        doc = yaml.safe_load(open(path))

        def walk(node, parent=None):
            if isinstance(node, dict):
                if node.get('name') == task_name and isinstance(node.get('attributes'), list):
                    for a in node['attributes']:
                        if a.get('name') == 'temperature' and 'default' in a:
                            return float(a['default'])
                for v in node.values():
                    r = walk(v, node)
                    if r is not None:
                        return r
            elif isinstance(node, list):
                for v in node:
                    r = walk(v, parent)
                    if r is not None:
                        return r
            return None
        return walk(doc)

    def test_membrane_equilibrate_default_still_310(self):
        found = self._default_for('membrane_equilibrate')
        self.assertIsNotNone(found, 'POSITIVE CONTROL: no membrane_equilibrate temperature default '
                                    'found in the schema, so this test proves nothing')
        self.assertEqual(found, _MEMBRANE_EQUILIBRATE_DEFAULT_TEMPERATURE)

    def test_the_diagnostic_exposes_a_temperature_override(self):
        """A user must be able to say what they mean; before the fix there was no such key."""
        import os
        import pestifer
        doc = yaml.safe_load(open(os.path.join(os.path.dirname(pestifer.__file__),
                                               'schema', 'base.yaml')))

        def find(node):
            if isinstance(node, dict):
                if node.get('name') == 'diagnose_differential_stress':
                    return node
                for v in node.values():
                    r = find(v)
                    if r:
                        return r
            elif isinstance(node, list):
                for v in node:
                    r = find(v)
                    if r:
                        return r
        block = find(doc)
        self.assertIsNotNone(block)
        names = [a['name'] for a in block['attributes']]
        self.assertIn('temperature', names)


class TestTheDiagnosticUsesIt(unittest.TestCase):
    """Aimed at the CALL SITE.  A correct helper that nobody calls is the original bug."""

    def setUp(self):
        self.src = inspect.getsource(MakeMembraneSystemTask.diagnose_differential_stress)

    def test_the_protocol_carries_a_temperature(self):
        line = next((ln for ln in self.src.splitlines()
                     if 'ensemble=' in ln and 'NPT' in ln), None)
        self.assertIsNotNone(line, 'POSITIVE CONTROL: no NPT stage found in the pass')
        self.assertIn('temperature=temperature', line,
                      f'the sampling stage does not set a temperature: {line.strip()!r}')

    def test_the_temperature_is_inherited_not_constant(self):
        self.assertIn('_relaxation_temperature(', self.src,
                      'the pass no longer reads the temperature the membrane equilibrated at')
        self.assertIn("get('quilt'", self.src,
                      'the inherited temperature is not taken from the QUILT protocol')

    def test_an_explicit_setting_overrides_the_inherited_one(self):
        self.assertIn("diag_specs.get('temperature')", self.src)


class TestTheDiagnosticDoesNotAdvanceTheMembrane(unittest.TestCase):

    def test_the_pass_is_run_read_only(self):
        src = inspect.getsource(MakeMembraneSystemTask.diagnose_differential_stress)
        call = src[src.index('self.equilibrate_bilayer('):]
        self.assertIn('preserve_state=True', call[:400],
                      'the diagnostic pass still rekeys its own final frame over quilt_state')

    def test_preserve_state_returns_before_the_rekey(self):
        """Aimed at the mechanism, not the flag.

        The flag could exist and be read after the damage was done.  The guard has to sit
        ahead of the rekey, which is the statement that replaces the embedded structure.
        """
        src = inspect.getsource(MakeMembraneSystemTask.equilibrate_bilayer)
        self.assertIn('preserve_state', src)
        guard = src.index('if preserve_state')
        rekey = src.index("rekey('state'")
        self.assertLess(guard, rekey, 'the preserve_state guard is after the rekey it must prevent')
        ret = src.index('return', guard)
        self.assertLess(ret, rekey, 'preserve_state does not return before the rekey')

    def test_the_other_callers_are_untouched(self):
        """What must NOT change.

        The patch and quilt relaxations advance the membrane -- that is their job -- so only the
        diagnostic may pass the flag.  A helper whose new behaviour leaked into its existing
        callers would be the worse bug.
        """
        src = inspect.getsource(MakeMembraneSystemTask)
        calls = [ln for ln in src.splitlines() if 'self.equilibrate_bilayer(' in ln]
        # patch, gridded quilt, tiled quilt, and the diagnostic
        self.assertEqual(len(calls), 4, f'expected 4 call sites, found {len(calls)}')
        self.assertEqual(sum(1 for ln in calls if 'preserve_state' in ln), 0,
                         'preserve_state should appear on the continuation line of the '
                         'diagnostic call only')
        self.assertEqual(src.count('preserve_state=True'), 1,
                         'exactly one caller may run the helper read-only')

    def test_default_is_still_advancing(self):
        sig = inspect.signature(MakeMembraneSystemTask.equilibrate_bilayer)
        self.assertIs(sig.parameters['preserve_state'].default, False)


if __name__ == '__main__':
    unittest.main()


class TestTheUnconfiguredCaseIsDeterminate(unittest.TestCase):
    """An absent quilt protocol is not an unknown, and treating it as one was my first fix's bug.

    `equilibrate_bilayer` substitutes `default_relaxation_protocol()` when handed none, so an
    unconfigured quilt relaxes under all-`md` stages at 300 K.  My first version warned and
    guessed 310 K there -- the `membrane_equilibrate` default -- which would have reintroduced
    the very 10 K mismatch this change exists to remove, in the one configuration where nobody
    had written a temperature to compare against.  An existing test caught it.
    """

    def test_the_hard_coded_fallback_runs_at_the_md_default(self):
        from pestifer.tasks.make_membrane_system import default_relaxation_protocol
        self.assertEqual(_T(default_relaxation_protocol()), _MD_DEFAULT_TEMPERATURE)

    def test_the_fallback_is_a_fresh_list_each_call(self):
        """NEGATIVE CONTROL: the stage loop writes into these dicts in place.

        A shared module-level list would carry one bilayer's `addl_paramfiles` and
        `other_parameters` into the next one's protocol.
        """
        from pestifer.tasks.make_membrane_system import default_relaxation_protocol
        a, b = default_relaxation_protocol(), default_relaxation_protocol()
        self.assertIsNot(a, b)
        self.assertIsNot(a[0]['md'], b[0]['md'])
        a[0]['md']['addl_paramfiles'] = ['leaked.prm']
        self.assertNotIn('addl_paramfiles', default_relaxation_protocol()[0]['md'])

    def test_one_definition_serves_both_readers(self):
        """The runner and the temperature lookup must not hold separate copies of the protocol."""
        src = inspect.getsource(MakeMembraneSystemTask.equilibrate_bilayer)
        self.assertIn('default_relaxation_protocol()', src)
        diag = inspect.getsource(MakeMembraneSystemTask.diagnose_differential_stress)
        self.assertIn('default_relaxation_protocol()', diag)
        self.assertNotIn("{'md': dict(ensemble='NPAT'", src,
                         'equilibrate_bilayer still carries an inline copy of the fallback')
