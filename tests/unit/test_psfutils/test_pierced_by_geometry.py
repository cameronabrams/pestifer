# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""`PSFRing.pierced_by` on synthetic geometry: the bond SEGMENT must cross the ring plane INSIDE
the ring, wherever inside.

Until 3.29.5 the test gated on distance: a bond whose midpoint lay more than half a bond length
from the ring COM was rejected, under a reason string claiming "both atoms on same side of ring
plane".  So every off-center threading was invisible.  It went unnoticed while threadings came
out of the build roughly centered (and with bonds stretched long by the minimizer, which widened
the gate), and it bit as soon as the resolver began making small moves: on example 17's
replica-1 quilt the 3.29.4 resolver reported 6 of 6 piercings cleared, and 3 were still threaded,
crossing 1-2 A off-center.  The resolver judges its own poses with this function, so it
certified its own non-fixes; the terminate gate shared the blind spot.

Pure geometry, no files, so it runs on CI.
"""
import unittest

import numpy as np

from pestifer.psfutil.psfring import PSFRing
from pestifer.psfutil.psfbond import PSFBond

CC = 1.40          # aromatic C-C, so the hexagon's circumradius is 1.40 A
HALF = 0.765       # half a 1.53 A C-C bond -- the width of the gate that used to apply


def ring(radius=CC, n=6, tilt=0.0):
    t = np.linspace(0, 2 * np.pi, n, endpoint=False)
    P = np.stack([radius * np.cos(t), radius * np.sin(t), tilt * np.sin(2 * t)], axis=1)
    r = PSFRing(list(range(1, n + 1)))
    r.P = P
    r.calculate_stuff()
    return r


def bond(a, b):
    x = PSFBond([101, 102])
    x.P = np.array([a, b], dtype=float)
    x.calculate_stuff()
    return x


def through(x, y, half=HALF):
    """A bond along z, crossing the ring plane at (x, y)."""
    return bond([x, y, -half], [x, y, +half])


class TestPiercedBy(unittest.TestCase):

    def test_a_centered_threading_is_pierced(self):
        self.assertTrue(ring().pierced_by(through(0.0, 0.0))['pierced'])

    def test_an_OFF_CENTER_threading_is_pierced(self):
        """The regression.  1.0 A off-center is past the old half-bond gate (0.765 A) but well
        inside the hexagon (inscribed radius 1.21 A).  The old code returned not-pierced."""
        r, b = ring(), through(1.0, 0.0)
        self.assertGreater(np.linalg.norm(b.COM - r.COM), HALF,
                           'this case must sit beyond the old gate, or it tests nothing')
        self.assertTrue(r.pierced_by(b)['pierced'])

    def test_off_center_along_a_vertex_direction_too(self):
        self.assertTrue(ring().pierced_by(through(0.6, 0.9))['pierced'])

    def test_a_bond_crossing_the_plane_outside_the_ring_is_not_pierced(self):
        res = ring().pierced_by(through(2.0, 0.0))
        self.assertFalse(res['pierced'])

    def test_a_bond_wholly_on_one_side_is_not_pierced(self):
        """Pointing straight at the hole from above -- the LINE threads the ring, the SEGMENT
        does not.  Removing the gate without adding a real straddle test would get this wrong."""
        res = ring().pierced_by(bond([0, 0, 0.3], [0, 0, 1.83]))
        self.assertFalse(res['pierced'])
        self.assertEqual(res['reason'], 'both atoms on same side of ring plane')

    def test_a_bond_lying_flat_above_the_ring_is_not_pierced(self):
        self.assertFalse(ring().pierced_by(bond([-0.7, 0, 1.0], [0.7, 0, 1.0]))['pierced'])

    def test_a_puckered_ring_still_works(self):
        """Sugar chairs and sterol A/B rings are not planar; the straddle test uses a mean
        plane."""
        r = ring(tilt=0.25)
        self.assertTrue(r.pierced_by(through(0.9, 0.0))['pierced'])
        self.assertFalse(r.pierced_by(through(2.2, 0.0))['pierced'])

    def test_a_five_membered_ring(self):
        r = ring(radius=1.19, n=5)                       # cyclopentane-like, as in a sterol D-ring
        self.assertTrue(r.pierced_by(through(0.8, 0.0))['pierced'])
        self.assertFalse(r.pierced_by(through(1.6, 0.0))['pierced'])


if __name__ == '__main__':
    unittest.main()
