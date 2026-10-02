# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
The axial penalty that keeps acyl tails on the tail side of the headgroup.

Without it the MC sampler's confinement is a cylinder about the membrane normal -- infinite along
that normal -- so a chain can fold 180 degrees back past the headgroup and out of the bilayer.
The trans-ordering field cannot object: its per-C-H term is ``1/2(3cos^2(theta) - 1)``, which
depends on ``cos^2`` and so is invariant under ``theta -> 180 - theta``.  An inverted chain scores
EXACTLY as well as an extended one, which made folding a free way to satisfy the bias.  Measured
2026-10-02: 10/10 POPE and 8/10 SOPE/SOPS drawn conformers had a tail tip above their own
phosphate, and the leaflets built 6.6-13 A too thin as a result.

It is a penalty and not a hard ceiling on purpose, and that is what these tests pin.
"""
import numpy as np
import unittest

from pestifer.charmmff.athermal_mc import MoleculeMC, RotatableBond, run_mc


def _two_bead_chain(z_top=5.0):
    """A minimal molecule: one rotatable bond swinging a 'tail' bead about a pivot."""
    coords = np.array([[0.0, 0.0, 10.0],    # 0: head reference
                       [0.0, 0.0, 0.0],     # 1: pivot anchor
                       [1.5, 0.0, 0.0],     # 2: pivot tip
                       [3.0, 0.0, z_top]])  # 3: the moving 'tail' bead
    radii = np.full(4, 0.1)
    bond = RotatableBond(a=1, b=2, moving=np.array([False, False, False, True]), i=0, j=3)
    confined = np.array([False, False, False, True])
    return coords, radii, bond, confined


class TestAxialPenaltyIsNotAWall(unittest.TestCase):

    def _mol(self, limit, z_top=5.0):
        coords, radii, bond, confined = _two_bead_chain(z_top)
        return MoleculeMC(coords=coords, radii=radii, rotatable=[bond],
                          axis_point=np.zeros(3), axis_dir=np.array([0.0, 0.0, 1.0]),
                          cylinder_radius=float('inf'), confined=confined,
                          axial_limit=limit)

    def test_zero_bias_leaves_the_sampler_untouched(self):
        """The default must reproduce the old sampler exactly, bit for bit."""
        a = run_mc(self._mol(0.0), nsamples=4, n_equil=50, n_decorr=10, seed=3, axial_bias=0.0)
        b = run_mc(self._mol(float('inf')), nsamples=4, n_equil=50, n_decorr=10, seed=3)
        for x, y in zip(a, b):
            np.testing.assert_allclose(x, y)

    def test_the_penalty_pushes_the_tail_below_the_limit(self):
        biased = run_mc(self._mol(0.0), nsamples=8, n_equil=400, n_decorr=20, seed=5,
                        axial_bias=5.0)
        free = run_mc(self._mol(0.0), nsamples=8, n_equil=400, n_decorr=20, seed=5,
                      axial_bias=0.0)
        mean_biased = np.mean([s[3, 2] for s in biased])
        mean_free = np.mean([s[3, 2] for s in free])
        self.assertLess(mean_biased, mean_free,
                        'the penalty did not lower the confined atom')

    def test_uphill_moves_are_accepted_sometimes(self):
        """THE defining property, and the one a hard wall does not have.

        Start the tail BELOW the limit and apply a weak penalty.  A penalty makes an excursion
        above the limit expensive but reachable, so some sample must exceed it.  A wall -- of any
        flavour, including the monotonically-inward one tried first -- makes it unreachable from
        inside, so no sample ever can.

        This replaces an earlier test that asserted a violating start could descend.  Mutation
        testing killed it: turning the penalty back into a wall left all five tests green, because
        a two-bead toy can always descend monotonically and so never reproduces the trap that hurt
        SOPE on real lipids.  A control that cannot fail is the bug it is guarding against.
        """
        mol = self._mol(0.0, z_top=-5.0)         # starts BELOW the limit
        samples = run_mc(mol, nsamples=40, n_equil=200, n_decorr=10, seed=17, axial_bias=0.05)
        above = max(float(s[3, 2]) for s in samples)
        self.assertGreater(above, 0.0,
                           'no sample ever rose above the limit -- the penalty is acting as a '
                           'wall, which traps a conformer that starts folded')

    def test_a_violating_start_can_still_escape(self):
        """THE POINT.  A hard ceiling traps a conformer that starts above the limit, because the
        path to the unfolded basin runs up and over it.

        Measured on real lipids 2026-10-02: under a hard ceiling SOPE went from 8/10 folded to
        10/10, and DMPC -- which never folds -- lost 2.06 A of extension to the blocked pivots.
        A penalty leaves every move reachable, just expensive, so the walk can leave the basin it
        started in.  If this ever becomes a wall again, this test fails.
        """
        mol = self._mol(0.0, z_top=5.0)          # starts ABOVE the limit
        samples = run_mc(mol, nsamples=10, n_equil=600, n_decorr=20, seed=11, axial_bias=2.0)
        below = [s for s in samples if s[3, 2] <= 0.0]
        self.assertTrue(below, 'no sample ever got below the limit -- the penalty acts as a wall')

    def test_an_infinite_limit_disables_it_regardless_of_bias(self):
        a = run_mc(self._mol(float('inf')), nsamples=4, n_equil=50, n_decorr=10, seed=7,
                   axial_bias=9.0)
        b = run_mc(self._mol(float('inf')), nsamples=4, n_equil=50, n_decorr=10, seed=7)
        for x, y in zip(a, b):
            np.testing.assert_allclose(x, y)


class TestBuildLipidMCDerivesTheLimit(unittest.TestCase):

    def test_the_limit_comes_from_the_lowest_head_reference_atom(self):
        """For a glycerolipid `heads` is [N, O22, O32] -- the headgroup nitrogen and the two ester
        oxygens.  The minimum is the ester plane, where the chains actually attach; a mean could be
        lifted by a leaning choline."""
        from pestifer.charmmff.athermal_mc import build_lipid_mc
        coords = np.array([[0.0, 0.0, 12.0],    # head N (high)
                           [0.0, 0.0, 4.0],     # ester O (the real ceiling)
                           [0.0, 0.0, 2.0],
                           [1.5, 0.0, 1.0],
                           [3.0, 0.0, 0.0],
                           [4.5, 0.0, -1.0]])
        masses = [14.0, 16.0, 12.0, 12.0, 12.0, 12.0]
        bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)]
        mol = build_lipid_mc(coords, None, masses, bonds, head_indices=[0, 1],
                             tail_indices=[5], rmin_half=[1.0] * 6, cylinder_radius=10.0)
        axial = np.dot(coords[[0, 1]] - mol.axis_point, mol.axis_dir)
        self.assertAlmostEqual(mol.axial_limit, float(axial.min()), places=9)
        self.assertLess(mol.axial_limit, float(axial.max()),
                        'POSITIVE CONTROL: the two head atoms must differ, or this proves nothing')


if __name__ == '__main__':
    unittest.main()
