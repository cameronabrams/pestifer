import unittest

import numpy as np

from pestifer.molecule.graft_geom import guess_bonds, pendant_indices


class TestGuessBonds(unittest.TestCase):
    def test_linear_chain(self):
        # four carbons ~1.5 A apart along x -> a bonded chain 0-1-2-3
        coords = np.array([[0.0, 0, 0], [1.5, 0, 0], [3.0, 0, 0], [4.5, 0, 0]])
        adj = guess_bonds(coords, ['C', 'C', 'C', 'C'])
        self.assertEqual(adj[0], {1})
        self.assertEqual(adj[1], {0, 2})
        self.assertEqual(adj[2], {1, 3})
        self.assertEqual(adj[3], {2})

    def test_non_bonded_gap_not_bonded(self):
        # a 3 A gap is well beyond a C-C bond -> no bond
        coords = np.array([[0.0, 0, 0], [3.0, 0, 0]])
        self.assertEqual(guess_bonds(coords, ['C', 'C']), [set(), set()])

    def test_c_o_and_c_h_bonds(self):
        # C at origin, O at 1.43 (bond), H at 1.09 (bond), a far O at 2.6 (no bond)
        coords = np.array([[0.0, 0, 0], [1.43, 0, 0], [0.0, 1.09, 0], [0.0, 0, 2.6]])
        adj = guess_bonds(coords, ['C', 'O', 'H', 'O'])
        self.assertEqual(adj[0], {1, 2})     # C-O and C-H
        self.assertNotIn(3, adj[0])          # the 2.6 A oxygen is not bonded


class TestPendantIndices(unittest.TestCase):
    def _tree_adj(self):
        # 0-1-2-3 with a branch 1-4-5 ; as adjacency sets
        adj = [set() for _ in range(6)]
        for a, b in [(0, 1), (1, 2), (2, 3), (1, 4), (4, 5)]:
            adj[a].add(b)
            adj[b].add(a)
        return adj

    def test_pendant_is_downstream_subtree(self):
        adj = self._tree_adj()
        # cut 1->2 : downstream of 2 is {3} (2 itself is index_j, excluded)
        self.assertEqual(pendant_indices(adj, 1, 2), {3})
        # cut 0->1 : downstream of 1 is everything except 0 and 1 -> {2,3,4,5}
        self.assertEqual(pendant_indices(adj, 0, 1), {2, 3, 4, 5})
        # cut 1->4 : downstream of 4 is {5}
        self.assertEqual(pendant_indices(adj, 1, 4), {5})

    def test_not_bonded_returns_empty(self):
        adj = self._tree_adj()
        self.assertEqual(pendant_indices(adj, 0, 3), set())   # 0 and 3 are not bonded

    def test_ring_does_not_infinite_loop(self):
        # a 3-cycle 0-1-2-0 plus pendant 2-3
        adj = [set() for _ in range(4)]
        for a, b in [(0, 1), (1, 2), (2, 0), (2, 3)]:
            adj[a].add(b)
            adj[b].add(a)
        # cutting 0->1, the BFS from 1 can reach 2 (via the ring) and 3, never crossing 0 or 1
        self.assertEqual(pendant_indices(adj, 0, 1), {2, 3})


if __name__ == '__main__':
    unittest.main()
