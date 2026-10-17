# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The rotation resolver's clash scorer must see across periodic boundaries.

`RingChecker.clash_count` scores a trial rotation by how many new heavy-atom contacts it makes.
It built its KD-tree without a `boxsize`, so a rotation that swung a tail into a neighbour
ACROSS a cell face scored zero while putting two heavy atoms on top of each other.

Measured 2026-10-16 on an example-17 patchB: `ring_check` cleared eight CHL1-ring piercings by
rotation, reported "0 new heavy-atom clash(es)" for each, and left `SOPE 18 C29` and
`SOPE 98 C26` **0.451 A** apart under the minimum image -- 65.7 A apart raw, which is what the
scorer saw.  The structure was then handed to dynamics.

    after solvate (pre-minimize)   14.148 A
    after the opening minimize     13.944 A
    after the ring_check rotations  0.451 A     <- created by a rotation scored as clean

This is the same blindness `_check_fast` had before d89cc518.  That fix went to the piercing
detector and not to the clash scorer beside it.
"""
import unittest
from types import SimpleNamespace

import numpy as np

from pestifer.psfutil.psfring import RingChecker


def _probe(n_atoms, neighbors=None):
    """Minimal stand-in carrying only what `clash_count` touches."""
    neighbors = neighbors or {}
    graph = SimpleNamespace(neighbors=lambda s: iter(neighbors.get(s, ())))
    return SimpleNamespace(_heavy=np.ones(n_atoms, dtype=bool),
                           _serials=list(range(1, n_atoms + 1)),
                           _row_of_serial={s: s - 1 for s in range(1, n_atoms + 1)},
                           topol=SimpleNamespace(G=graph))


class TestClashCountAcrossPeriodicBoundaries(unittest.TestCase):

    # The 3x3 cell, which is what every caller actually holds -- `check_coords(scratch, box,
    # ...)` sits on the line above each `clash_count` call and takes the matrix.  The first
    # version of this test passed a bare 3-vector, so it exercised a shape no caller uses and
    # went green while the real call raised
    # `ValueError: operands could not be broadcast together with shapes (15189,3) (3,3)`
    # on the first patchB rotation of a real build.
    BOX = np.diag([20.0, 20.0, 20.0])
    EDGES = np.array([20.0, 20.0, 20.0])

    def _coords(self, gap):
        """Two atoms either side of the x face, `gap` apart under the minimum image."""
        return np.array([[0.2, 10.0, 10.0], [0.2 - gap + self.EDGES[0], 10.0, 10.0]])

    def test_a_contact_across_a_face_is_counted_when_the_box_is_given(self):
        """The defect in one assertion: this used to be 0."""
        n = RingChecker.clash_count(_probe(2), self._coords(0.451), moved_rows=[0],
                                    cutoff=2.0, box=self.BOX)
        self.assertEqual(n, 1, 'a 0.451 A contact across a cell face was not counted')

    def test_the_same_pair_is_invisible_without_a_box(self):
        """NEGATIVE CONTROL for the test itself: it must be the box that makes the difference,
        not the geometry.  Raw separation here is 19.5 A."""
        n = RingChecker.clash_count(_probe(2), self._coords(0.451), moved_rows=[0], cutoff=2.0)
        self.assertEqual(n, 0)

    def test_a_genuinely_distant_pair_is_not_counted_either_way(self):
        coords = np.array([[2.0, 10.0, 10.0], [12.0, 10.0, 10.0]])
        for box in (self.BOX, self.EDGES, None):
            with self.subTest(box=box):
                self.assertEqual(
                    RingChecker.clash_count(_probe(2), coords, moved_rows=[0],
                                            cutoff=2.0, box=box), 0)

    def test_an_ordinary_in_cell_contact_is_still_counted(self):
        coords = np.array([[10.0, 10.0, 10.0], [10.6, 10.0, 10.0]])
        self.assertEqual(RingChecker.clash_count(_probe(2), coords, moved_rows=[0],
                                                 cutoff=2.0, box=self.BOX), 1)

    def test_a_bonded_neighbor_is_still_exempt(self):
        """The hinge partner is an expected close contact, periodic or not."""
        coords = np.array([[10.0, 10.0, 10.0], [10.6, 10.0, 10.0]])
        probe = _probe(2, neighbors={1: (2,), 2: (1,)})
        self.assertEqual(RingChecker.clash_count(probe, coords, moved_rows=[0],
                                                 cutoff=2.0, box=self.BOX), 0)


    def test_a_bare_edge_vector_works_too(self):
        n = RingChecker.clash_count(_probe(2), self._coords(0.451), moved_rows=[0],
                                    cutoff=2.0, box=self.EDGES)
        self.assertEqual(n, 1)

    def test_the_cell_matrix_is_the_shape_callers_pass(self):
        """Guards the regression directly: a 3x3 must not be broadcast against (N,3)."""
        self.assertEqual(np.asarray(self.BOX).shape, (3, 3))
        n = RingChecker.clash_count(_probe(2), self._coords(0.451), moved_rows=[0],
                                    cutoff=2.0, box=self.BOX)
        self.assertEqual(n, 1)

    def test_every_call_site_passes_a_box(self):
        """A periodic-aware scorer helps nobody if a caller leaves the box out -- the same
        shape as the cell_from_xsc sweep that cleared a wrapper and missed the layer below."""
        import inspect, re
        from pestifer.psfutil import ring_resolve
        from pestifer.tasks import ringcheck
        for mod in (ring_resolve, ringcheck):
            for call in re.findall(r'clash_count\(([^)]*)\)', inspect.getsource(mod)):
                self.assertIn('box=', call, f'{mod.__name__}: clash_count({call}) has no box')
