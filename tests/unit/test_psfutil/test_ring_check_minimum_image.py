# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""The fast targeted ring check sees a piercing across a periodic boundary.

``RingChecker._check_fast`` narrows the bonds it tests by ``|bond_midpoint - ring_COM| <=
cutoff``. That distance was taken on **raw** coordinates, and ``mic_shift`` was applied only to
bonds that had already survived the filter -- so a bond lying across a periodic boundary from the
ring measured a box length away and was discarded before it could be shifted.

Measured on an ex17 calibration patch (2026-10-13, box 76 x 66 x 86 A): a real piercing of
``CHL1 L1-199`` sat **1.13 A** from the ring under the minimum image and **66.56 A** in raw
coordinates. ``_check_fast`` reported the ring unpierced; the image-aware full ``check`` reported
it pierced.

The docstring asserted this path was "exact, not approximate". It was exact only for a piercing
wholly inside the box -- and in a membrane, lipids straddle the cell edge by construction, so the
boundary case is the common one.

**Why it mattered far beyond a missed report.** The rotation resolver's acceptance test is
``check_coords(candidate, box, [target])``. For a boundary-crossing piercing that returned "not
pierced" whatever the candidate did, so the FIRST rotation tried was accepted and logged as
cleared while the piercing survived -- to be found by the next full check and resolved by
deleting a cholesterol, whose threading partner was then stranded in a contorted pose and
RATTLE-failed at the first dynamics step. One observed build death traces along exactly that
chain (job 26469605 rep-05, atom 30728 = POPC L1-292).

Verified end to end on the real structures: with this fix the quilt stage resolves every piercing
by rotation with no deletions, where production deleted 4-5 cholesterols per build.
"""
import unittest

import numpy as np

from pestifer.psfutil.psfring import RingChecker


class _Ring:
    """Minimal stand-in for a PSFRing: a planar hexagon at a chosen center."""

    def __init__(self, center, rows):
        self.COM = np.asarray(center, dtype=float)
        self._rows = rows
        self.idx_list = [100 + i for i in range(len(rows))]
        self.segname, self.resid = 'L1', 199
        self.P = None
        self.pierced = []

    def calculate_stuff(self):
        pass

    def pierced_by(self, bond):
        # every bond that reaches this method counts as a piercing, so the test isolates the
        # CANDIDATE FILTER -- which is where the defect was -- rather than the geometry test
        self.pierced.append(bond)
        return {'pierced': True}


class TestTheCandidateFilterUsesTheMinimumImage(unittest.TestCase):

    def _checker(self, bond_mid, ring_com, box_edges):
        ck = RingChecker.__new__(RingChecker)
        ck.cutoff = 4.0
        half = np.array([0.4, 0.0, 0.0])
        coords = np.array([bond_mid - half, bond_mid + half,
                           *[ring_com] * 6], dtype=float)
        ck._bond_rows = np.array([[0, 1]])
        ring = _Ring(ring_com, np.arange(2, 8))

        class _B:
            idx_list = [1, 2]
            P = None

            def calculate_stuff(self):
                pass

            def mic_shift(self, com, box):
                return self
        ck._bonds = [_B()]
        ck._ring_index = {('L1', '199'): [ring]}
        ck._piercespec = lambda r, b, i: {'piercee': {'segname': 'L1', 'resid': 199},
                                          'piercer': {'segname': 'L1', 'resid': 4}}
        box = np.diag(np.asarray(box_edges, dtype=float))
        return ck, ring, coords, box

    def test_a_bond_across_the_boundary_IS_found(self):
        """The regression. Raw separation 66.5 A, minimum image 1.1 A, cutoff 4.0."""
        edges = [75.98, 65.80, 86.22]
        ck, ring, coords, box = self._checker(bond_mid=[1.0, 1.0, 1.0],
                                              ring_com=[74.9, 64.9, 85.1], box_edges=edges)
        raw = np.linalg.norm(np.array([1.0, 1.0, 1.0]) - np.array([74.9, 64.9, 85.1]))
        self.assertGreater(raw, 4.0, 'POSITIVE CONTROL: the raw separation must exceed the cutoff')
        out = ck._check_fast(coords, box, [('L1', 199)])
        self.assertEqual(len(out), 1, 'the boundary-crossing bond was filtered out again')

    def test_a_bond_genuinely_far_away_is_still_rejected(self):
        """NEGATIVE CONTROL: the fix must not turn the filter into "accept everything"."""
        edges = [200.0, 200.0, 200.0]
        ck, ring, coords, box = self._checker(bond_mid=[10.0, 10.0, 10.0],
                                              ring_com=[90.0, 90.0, 90.0], box_edges=edges)
        out = ck._check_fast(coords, box, [('L1', 199)])
        self.assertEqual(out, [], 'a bond 138 A away (even under the minimum image) was accepted')

    def test_a_bond_close_in_raw_coordinates_is_still_found(self):
        """The case that always worked must keep working."""
        edges = [200.0, 200.0, 200.0]
        ck, ring, coords, box = self._checker(bond_mid=[50.0, 50.0, 50.0],
                                              ring_com=[51.0, 50.0, 50.0], box_edges=edges)
        out = ck._check_fast(coords, box, [('L1', 199)])
        self.assertEqual(len(out), 1)

    def test_a_non_periodic_run_is_unaffected(self):
        """box=None must not raise and must behave as the raw-distance filter."""
        ck, ring, coords, _box = self._checker(bond_mid=[50.0, 50.0, 50.0],
                                               ring_com=[51.0, 50.0, 50.0], box_edges=[1, 1, 1])
        self.assertEqual(len(ck._check_fast(coords, None, [('L1', 199)])), 1)
        ck2, _r, coords2, _b = self._checker(bond_mid=[0.0, 0.0, 0.0],
                                             ring_com=[80.0, 0.0, 0.0], box_edges=[1, 1, 1])
        self.assertEqual(ck2._check_fast(coords2, None, [('L1', 199)]), [],
                         'without a cell there is no image to fold into')


if __name__ == '__main__':
    unittest.main()
