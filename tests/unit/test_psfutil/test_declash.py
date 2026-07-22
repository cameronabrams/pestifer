import unittest

import numpy as np
from scipy.spatial import cKDTree

from pestifer.psfutil.declash import declash_pendant, declash_loop, _contact_count


class TestDeclashPendant(unittest.TestCase):
    def test_clears_a_resolvable_clash(self):
        # mover atom sits on top of an env atom; a rotation about the i-j bond (the z axis
        # through the origin) swings it away
        coords = np.array([[0., 0, 0], [0, 0, 3], [1, 0, 3], [1, 0, 3]])   # i, j, mover, env
        rng = np.random.default_rng(0)
        n = declash_pendant(coords, np.array([2]), np.array([3]),
                            [(0, 1)], [np.array([2])], 50, 1.5, rng)
        self.assertEqual(n, 0)
        self.assertGreater(np.linalg.norm(coords[2] - coords[3]), 1.5)

    def test_greedy_never_increases_clashes(self):
        rng = np.random.default_rng(1)
        coords = rng.normal(scale=1.5, size=(20, 3))
        pend, env = np.arange(10), np.arange(10, 20)
        before = _contact_count(cKDTree(coords[env]), coords[pend], 2.0)
        declash_pendant(coords, pend, env, [(0, 1)], [np.array([2, 3, 4, 5])], 40, 2.0, rng)
        after = _contact_count(cKDTree(coords[env]), coords[pend], 2.0)   # env is static
        self.assertLessEqual(after, before)

    def test_no_bonds_is_a_noop(self):
        coords = np.random.default_rng(2).normal(size=(6, 3))
        orig = coords.copy()
        declash_pendant(coords, np.array([0, 1]), np.array([2, 3]), [], [], 10, 2.0,
                        np.random.default_rng(0))
        self.assertTrue(np.array_equal(coords, orig))


class TestDeclashLoop(unittest.TestCase):
    def test_reduces_loop_clash(self):
        # a "loop" mover atom overlaps an env atom; a phi rotation about N->CA moves it off
        coords = np.array([[0., 0, 0], [0, 0, 1], [1, 0, 2], [1, 0, 2]])   # N, CA, mover, env
        residues = [{'pivot': 0, 'axis_to': 1, 'movers': np.array([2]), 'frag_heavy': np.array([2])}]
        rng = np.random.default_rng(0)
        n = declash_loop(coords, residues, np.array([3]), 200, 1.5, rng)
        self.assertEqual(n, 0)
        self.assertGreater(np.linalg.norm(coords[2] - coords[3]), 1.5)

    def test_no_env_is_a_noop(self):
        coords = np.random.default_rng(3).normal(size=(5, 3))
        orig = coords.copy()
        residues = [{'pivot': 0, 'axis_to': 1, 'movers': np.array([2, 3]), 'frag_heavy': np.array([2, 3])}]
        declash_loop(coords, residues, np.array([], dtype=int), 10, 2.0, np.random.default_rng(0))
        self.assertTrue(np.array_equal(coords, orig))


if __name__ == '__main__':
    unittest.main()
