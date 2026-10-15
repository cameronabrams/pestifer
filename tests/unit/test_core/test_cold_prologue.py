# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""No hot, long-timestep, constrained dynamics on coordinates that were just built or moved.

A RATTLE constraint failure kills a build outright, and the way pestifer earns one is always
the same shape: a stage that modified coordinates is followed directly by dynamics that is hot
enough, long-stepped enough and constrained enough to turn residual strain into a divergent
step.  Two places in this repository already say so in prose:

- ``make_membrane_system._guard_pierced_lipids``: a threading "survives further minimization
  (it is topological, not just strained) and then RATTLE-fails the first dynamics step".  Its
  cure is ``ring_check`` THEN ``minimize``, spliced before the first non-minimize stage.
- example 25: "a full-temperature NVT start RATTLE-fails on the first step... a 50 K, 0.5 fs
  stage with rigidBonds off (no constraints to fail) lets the residual strain dissipate".

Example 17 had neither: its post-embed ladder minimized BEFORE ``ring_check`` and then went
straight to 410 K at 2 fs with ``rigidbonds all`` and freshly assigned 410 K velocities -- the
hottest stage in the example set, on the most complex system in it, on just-built coordinates.
These tests pin the repaired ladder, because a protocol is a claim like any other and nothing
else here would fail if it were reverted.

Pure config resolution: no VMD, no NAMD, no network.  Deliberately NOT marked ``needs_tools``
-- it has to run on CI, which is the only place that sees a fresh checkout.
"""
import unittest
from importlib.resources import files as pkg_files
from pathlib import Path

from pestifer.core.config import Config

# Locate the shipped examples through the PACKAGE, never a relative path.  `conftest.py`'s
# `change_test_dir` chdirs per test function, while `setUpClass` runs before the first such
# chdir -- so a relative path here resolves in `setUpClass` and raises inside a test body, and
# the module half-passes depending on where each load happens to sit.  Measured, not reasoned:
# the first version of this file did exactly that.
_EXAMPLES = Path(str(pkg_files('pestifer.resources'))) / 'examples'
_EX17 = str(_EXAMPLES / '17' / 'inputs' / 'hiv-mpertm3-membrane2.yaml')

# Above this, a stage is "hot" enough that residual strain matters.  310 K stages are routine
# and run everywhere; the melt is the outlier.
_HOT_K = 320.0


def _tasks(userfile):
    return Config(userfile=userfile, quiet=True).configure_new()['user']['tasks']


def _kind(task):
    (k, v), = task.items()
    return k, v


class TestExample17ColdPrologue(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.tasks = _tasks(_EX17)
        cls.kinds = [_kind(t) for t in cls.tasks]

    def test_a_minimize_follows_ring_check_before_anything_else(self):
        """``ring_check`` MOVES atoms -- it rotates a pendant branch or a side chain to unthread
        a piercing -- and whatever it moves is unminimized until something minimizes it."""
        idx = [i for i, (k, _) in enumerate(self.kinds) if k == 'ring_check']
        self.assertTrue(idx, 'example 17 no longer runs a top-level ring_check')
        for i in idx:
            nxt = self.kinds[i + 1]
            self.assertEqual(nxt[0], 'md', f'task after ring_check at {i} is {nxt[0]}, not md')
            self.assertEqual(str(nxt[1].get('ensemble', '')).casefold(), 'minimize',
                             f'task after ring_check at {i} is {nxt[1].get("ensemble")}, '
                             'not a minimize')

    def test_the_hot_melt_is_not_the_first_dynamics_after_the_build(self):
        """The 410 K pass must be preceded by dynamics that is cold, short-stepped and
        unconstrained -- the stage that cannot produce a RATTLE failure because there are no
        constraints in it."""
        dyn = [(i, v) for i, (k, v) in enumerate(self.kinds)
               if k == 'md' and str(v.get('ensemble', '')).casefold() != 'minimize']
        self.assertTrue(dyn, 'example 17 runs no dynamics at all')
        first_i, first = dyn[0]
        self.assertLess(float(first.get('temperature', 0)), _HOT_K,
                        f'first dynamics stage (task {first_i}) runs at '
                        f'{first.get("temperature")} K')
        op = first.get('other_parameters') or {}
        self.assertLessEqual(float(op.get('timestep', 2.0)), 1.0,
                             'first dynamics stage does not shorten the timestep')
        self.assertEqual(str(op.get('rigidbonds', 'all')).casefold(), 'none',
                         'first dynamics stage still constrains bonds, so RATTLE can fail in it')

    def test_every_hot_stage_is_preceded_by_a_cooler_one(self):
        """Temperature must be approached from below, not assigned outright."""
        dyn = [(i, v) for i, (k, v) in enumerate(self.kinds)
               if k == 'md' and str(v.get('ensemble', '')).casefold() != 'minimize']
        temps = [float(v.get('temperature', 0)) for _, v in dyn]
        hot = [n for n, t in enumerate(temps) if t >= _HOT_K]
        self.assertTrue(hot, 'no hot stage left in example 17 -- update this test with the reason')
        for n in hot:
            self.assertGreater(n, 0, 'a hot stage is the first dynamics in the build')
            self.assertLess(temps[n - 1], temps[n],
                            f'stage at {temps[n]} K is not approached from below '
                            f'(previous dynamics stage was {temps[n-1]} K)')


class TestOtherParametersOverrideRatherThanDuplicate(unittest.TestCase):
    """``other_parameters`` merges into a case-sensitive dict of NAMD directives, so a key
    spelled differently from the built-in one does not override it -- it adds a SECOND
    directive, and NAMD rejects the duplicate.  example 25 carries that warning in a comment;
    this fails if example 17's spellings drift off the built-ins."""

    def test_timestep_and_rigidbonds_replace_the_solvated_defaults(self):
        c = Config(userfile=_EX17, quiet=True).configure_new()
        solvated = c['user']['namd']['solvated']
        softened = False
        for i, t in enumerate(c['user']['tasks']):
            k, v = _kind(t)
            op = (v.get('other_parameters') or {}) if k == 'md' else {}
            if not ({'timestep', 'rigidbonds'} & set(op)):
                continue
            softened = True
            added = set(op) - set(solvated)
            self.assertFalse(added, f'task {i} introduces NAMD directive(s) {added} instead of '
                                    'overriding the built-in spelling; NAMD will see duplicates')
        self.assertTrue(softened, 'no stage softens timestep/rigidbonds any more')


class TestEveryShippedExampleStillResolves(unittest.TestCase):
    """A protocol edit that passes `yaml.safe_load` can still be rejected by the schema."""

    def test_all_example_configs_resolve(self):
        files = sorted(_EXAMPLES.glob('*/inputs/*.yaml'))
        self.assertGreater(len(files), 25, f'no example configs under {_EXAMPLES}')
        for f in files:
            with self.subTest(example=f):
                Config(userfile=str(f), quiet=True).configure_new()
