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
        lifted by a leaning choline.

        A per-atom limit (each atom bounded by its graph-nearest head) was tried 2026-10-02 to
        reach PMCL1, whose cardiolipin arms attach at very different heights.  It did not fix
        PMCL1 and it REGRESSED DSPE__Lo from 0/10 folded at +20.39 A to 10/10 at -0.56, because
        bounding a tail atom by a high head reference is looser than the ester plane.  One scalar,
        from the lowest head reference, plus the deterministic pre-pass, fixes five of the six.
        """
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


class TestUnfoldingPrePass(unittest.TestCase):
    """A folded START cannot be escaped by the penalty, so it is relieved deterministically first.

    Five shipped lipids are built from internal coordinates with a chain already 13-20 A above
    their own headgroup.  Unfolding needs a concerted swing that costs axial excess on the way,
    and the penalty that keeps good conformers good is exactly what forbids that detour -- raising
    `axial_bias` from 1.0 to 5.0 to 20.0 made folding MORE common, not less.
    """

    def test_a_folded_start_is_relieved_before_sampling(self):
        from pestifer.charmmff.athermal_mc import relieve_axial_excess
        coords, radii, bond, confined = _two_bead_chain(z_top=6.0)
        mol = MoleculeMC(coords=coords, radii=radii, rotatable=[bond],
                         axis_point=np.zeros(3), axis_dir=np.array([0.0, 0.0, 1.0]),
                         cylinder_radius=float('inf'), confined=confined,
                         axial_limit=0.0)
        out = relieve_axial_excess(mol, coords)
        self.assertLess(out[3, 2], coords[3, 2], 'the pre-pass did not lower the folded atom')

    def test_it_is_a_no_op_without_a_limit(self):
        from pestifer.charmmff.athermal_mc import relieve_axial_excess
        coords, radii, bond, confined = _two_bead_chain(z_top=6.0)
        mol = MoleculeMC(coords=coords, radii=radii, rotatable=[bond],
                         axis_point=np.zeros(3), axis_dir=np.array([0.0, 0.0, 1.0]),
                         cylinder_radius=float('inf'), confined=confined)
        np.testing.assert_allclose(relieve_axial_excess(mol, coords), coords)


class TestArmPivotsForMultiArmLipids(unittest.TestCase):
    """The unfolding pre-pass may rotate a whole phosphatidyl arm; the MC still may not.

    A bond is a tail torsion only when its tip-side fragment is all carbon, which keeps the
    headgroup rigid by design.  A cardiolipin's arm carries phosphorus and oxygen, so the bond
    joining it to the central glycerol is excluded and the arm cannot be reoriented at all --
    PMCL1 is built with its A-arm ester at +6.64 A while both its phosphates sit at -2.31 and
    +2.24, which no acyl torsion can fix and which `--refic-idx` 1, 2 and 3 all reproduce.

    These pivots are for the pre-pass ALONE.  It is a one-time search for a better starting point,
    and the IC-built arm orientation is an arbitrary starting choice rather than a physical
    equilibrium, so correcting it is fixing the input, not biasing the output.
    """

    @staticmethod
    def _cardiolipin_like():
        """Central carbon with TWO arms, each `C-O-P-O-C` then an all-carbon chain.

        The arm-joining bonds carry O and P on their tip side, so they are NOT acyl torsions.
        """
        coords, masses, bonds, names = [], [], [], []

        def add(z, m, nm):
            coords.append([0.0, 0.0, float(z)]); masses.append(m); names.append(nm)
            return len(coords) - 1

        centre = add(0.0, 12.011, 'C2')
        for arm, sign in (('A', 1.0), ('B', -1.0)):
            o1 = add(sign * 1.0, 15.999, f'O{arm}1')
            ph = add(sign * 2.0, 30.974, f'P{arm}')
            o2 = add(sign * 3.0, 15.999, f'O{arm}2')
            bonds += [(centre, o1), (o1, ph), (ph, o2)]
            prev = o2
            for k in range(2, 10):                      # an 8-carbon all-carbon chain
                c = add(sign * (3.0 + k), 12.011, f'C{arm}{k}')
                bonds.append((prev, c)); prev = c
        return np.array(coords), masses, bonds

    def test_an_arm_is_a_pre_pass_pivot_but_not_an_mc_torsion(self):
        from pestifer.charmmff.athermal_mc import build_lipid_mc
        coords, masses, bonds = self._cardiolipin_like()
        mol = build_lipid_mc(coords, None, masses, bonds, head_indices=[0],
                             tail_indices=[12, 24], rmin_half=[1.0] * len(masses),
                             cylinder_radius=50.0)

        def moved(bond):
            return frozenset(int(x) for x in np.flatnonzero(bond.moving))

        mc_sets = {moved(b) for b in mol.rotatable}
        extra_sets = {moved(b) for b in (mol.extra_pivots or [])}
        self.assertTrue(extra_sets, 'no arm pivot was derived for a two-armed lipid')
        self.assertFalse(mc_sets & extra_sets,
                         'an arm pivot leaked into the MC torsions; the sampled degrees of '
                         'freedom must be unchanged')
        # at least one arm pivot must carry a phosphorus -- i.e. move a whole arm, not a sub-chain
        phos = {i for i, m in enumerate(masses) if 30.5 < m < 31.5}
        self.assertTrue(any(s & phos for s in extra_sets),
                        'no pivot moves a phosphorus, so no whole arm can be reoriented')

    def test_the_pre_pass_does_not_widen_what_the_mc_samples(self):
        """The invariant, checked AFTER the pre-pass has run.

        Comparing `rotatable` against `extra_pivots` at construction time cannot see this: caught
        by mutation testing, where folding the arm pivots into `mol.rotatable` inside
        `relieve_axial_excess` left all the other tests green.  A control aimed one level off is
        the same bug it is guarding against, so this one runs the pre-pass first.
        """
        from pestifer.charmmff.athermal_mc import build_lipid_mc, relieve_axial_excess
        coords, masses, bonds = self._cardiolipin_like()
        mol = build_lipid_mc(coords, None, masses, bonds, head_indices=[0],
                             tail_indices=[12, 24], rmin_half=[1.0] * len(masses),
                             cylinder_radius=50.0)
        mol.axial_limit = -1.0                      # force real work out of the pre-pass
        before = [frozenset(int(x) for x in np.flatnonzero(b.moving)) for b in mol.rotatable]
        relieve_axial_excess(mol, mol.coords)
        after = [frozenset(int(x) for x in np.flatnonzero(b.moving)) for b in mol.rotatable]
        self.assertEqual(before, after,
                         'the pre-pass changed mol.rotatable, so the MC would now sample arm '
                         'rotations too -- these pivots are for the pre-pass alone')

    def test_an_ordinary_two_chain_lipid_gains_nothing(self):
        """POSITIVE CONTROL in reverse: a lipid whose chains are already reachable by acyl
        torsions must not acquire extra pivots, or this is changing every lipid rather than the
        multi-arm ones."""
        from pestifer.charmmff.athermal_mc import build_lipid_mc
        coords = np.array([[0.0, 0.0, 10.0 - k] for k in range(12)])
        masses = [14.007, 15.999] + [12.011] * 10
        bonds = [(k, k + 1) for k in range(11)]
        mol = build_lipid_mc(coords, None, masses, bonds, head_indices=[0, 1],
                             tail_indices=[11], rmin_half=[1.0] * 12, cylinder_radius=50.0)
        phos = []
        for b in (mol.extra_pivots or []):
            phos.append(int(np.count_nonzero(b.moving)))
        self.assertEqual(mol.extra_pivots, [],
                         f'a simple two-chain lipid gained {len(mol.extra_pivots or [])} arm '
                         f'pivots (moving {phos}); only multi-arm lipids should')


if __name__ == '__main__':
    unittest.main()
