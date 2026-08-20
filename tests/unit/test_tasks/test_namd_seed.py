# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the NAMD RNG seed a build hands each MD invocation.

NAMD's own default seed is drawn from the clock, so before this existed two runs of one config
produced different trajectories and neither could be reproduced.  These pin the two properties
that matter: a build is reproducible, and its individual NAMD runs are not correlated with one
another.
"""

import unittest
from types import SimpleNamespace

from pestifer.tasks.mdtask import MDTask


def _task(seed, controller_index=0, index=0, subtaskcount=0):
    t = MDTask.__new__(MDTask)
    t.namd_global_config = {'seed': seed}
    t.controller_index = controller_index
    t.index = index
    t.subtaskcount = subtaskcount
    return t


class TestNAMDSeed(unittest.TestCase):

    def test_zero_means_let_namd_decide(self):
        self.assertEqual(_task(0)._namd_seed(), 0)

    def test_absent_key_means_let_namd_decide(self):
        t = _task(0)
        t.namd_global_config = {}
        self.assertEqual(t._namd_seed(), 0)

    def test_same_base_and_position_is_reproducible(self):
        self.assertEqual(_task(27021972, 0, 3, 1)._namd_seed(),
                         _task(27021972, 0, 3, 1)._namd_seed())

    def test_different_base_gives_a_different_stream(self):
        # this is what makes a replica a replica
        self.assertNotEqual(_task(27021972, 0, 3, 0)._namd_seed(),
                            _task(27021973, 0, 3, 0)._namd_seed())

    def test_consecutive_replica_seeds_are_not_adjacent(self):
        """Replicas are seeded 27021972/3/4 by both sweep scripts, so consecutive base seeds are
        the normal case rather than an edge one.

        The earlier derivation XORed the base against a per-stream constant, which leaves the
        base's low bits untouched: a real triplicate came out with NAMD seeds 1112820144,
        1112820145 and 1112820146.  Adjacent seeds are not necessarily independent streams for a
        generator seeded by simple state initialization, and it is not a property worth having to
        defend in review.  Hashing the base together with the stream key avalanches it.
        """
        seeds = [_task(b, 0, 6, 0)._namd_seed() for b in (27021972, 27021973, 27021974)]
        gaps = [abs(b - a) for a, b in zip(seeds, seeds[1:])]
        self.assertGreater(min(gaps), 1_000_000,
                           f'consecutive replicas got near-adjacent seeds: {seeds}')

    def test_a_replica_set_stays_reproducible(self):
        """Separation must not have cost determinism: the same base still gives the same seed."""
        for b in (27021972, 27021973, 27021974):
            self.assertEqual(_task(b, 0, 6, 0)._namd_seed(), _task(b, 0, 6, 0)._namd_seed())

    def test_runs_within_one_build_are_also_well_separated(self):
        seeds = [_task(27021972, 0, i, 0)._namd_seed() for i in range(3, 9)]
        self.assertEqual(len(set(seeds)), len(seeds))
        pairs = [abs(a - b) for i, a in enumerate(seeds) for b in seeds[i + 1:]]
        self.assertGreater(min(pairs), 1_000_000, f'runs share a seed neighbourhood: {seeds}')

    def test_each_invocation_in_a_build_gets_its_own_stream(self):
        # reusing one seed for every NAMD run in a build would correlate them
        seeds = {_task(27021972, 0, i, 0)._namd_seed() for i in range(12)}
        self.assertEqual(len(seeds), 12)

    def test_subtasks_differ_too(self):
        seeds = {_task(27021972, 0, 5, k)._namd_seed() for k in range(8)}
        self.assertEqual(len(seeds), 8)

    def test_subcontrollers_differ_too(self):
        seeds = {_task(27021972, c, 0, 0)._namd_seed() for c in range(6)}
        self.assertEqual(len(seeds), 6)

    def test_seed_is_a_positive_int_in_namd_range(self):
        for base in (1, 27021972, 2 ** 31 - 2):
            for i in range(20):
                s = _task(base, 0, i, 0)._namd_seed()
                self.assertIsInstance(s, int)
                self.assertGreater(s, 0)
                self.assertLess(s, 2 ** 31)


class TestSeedOverride(unittest.TestCase):
    """``run --seed`` is defined as setting ``namd.seed``; it must land in the user config so it
    reaches the scripters, ``--complete-config``, and any taskless subconfiguration."""

    @classmethod
    def setUpClass(cls):
        from pestifer.core.config import Config
        cls.Config = Config
        cls.RM = Config(quiet=True).configure_new().RM

    def test_override_lands_in_user_config(self):
        c = self.Config(userdict={}, quiet=True, RM=self.RM, seed_override=98765).configure()
        self.assertEqual(c['user']['namd']['seed'], 98765)
        self.assertEqual(c.taskless_subconfig()['user']['namd']['seed'], 98765)

    def test_zero_override_is_honored_not_treated_as_unset(self):
        c = self.Config(userdict={}, quiet=True, RM=self.RM, seed_override=0).configure()
        self.assertEqual(c['user']['namd']['seed'], 0)

    def test_default_is_a_fixed_seed_so_builds_reproduce(self):
        c = self.Config(userdict={}, quiet=True, RM=self.RM).configure()
        self.assertTrue(c['user']['namd']['seed'])
