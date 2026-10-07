# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A conformer whose chains do not descend is rejected at draw time, in ANY entry.

The collapse guard added in 3.26.0 compares a tuned ensemble against its own unbiased one, so it
is blind to base-phase entries -- a base entry IS the unbiased ensemble and has no control to be
measured against.  That is where the remaining damage sat: 123 of 2380 conformers in 45 entries of
the shipped 3.26.0 collection, with fully SATURATED chains descending as little as 0.01 of their
contour, which no kink excuses.

The measure is per chain, against that chain's OWN attachment carbon, normalised by its all-trans
contour.  Per chain because the failure is one chain folding while its partner extends: a
molecule-wide head-to-tail number averages the two and reports neither, which is exactly how the
shipped ``POPC__Lo`` passed a "0 folded" audit while half its chains were hairpins.

Three things were checked before the threshold was set, and the checks are what the numbers in
these tests are:

- **Unsaturation does not force a low score.**  Median r by (chain length, sp2 count) runs
  0.81/0.75/0.73/0.75/0.72/0.60 from saturated to twelve sp2 carbons -- ~0.15, against a cut
  sitting 0.4 below the median.  An early hypothesis that polyunsaturated chains would be
  condemned was wrong, and so was a follow-up guess that PLPI33's low score was its linoleoyl
  chain: its unsaturated 18C chain scores 0.75 and its SATURATED 16C chain 0.46.
- **The cut is geometric, not tuned.**  Classifying every suspect entry's worst conformer by
  axial drop and radial displacement separates at 0.21/0.26: at or below 0.21 each one is a fold
  or a splay, at or above 0.26 each one descends.
- **Multi-arm lipids are not unfairly punished.**  Four-chain molecules are rejected more often
  (15% of conformers against 4%), but their PER-CHAIN rate is 3.7% against 2.2%, and 3.7% over
  four chains predicts 13.8% -- so it is order statistics, not a bias.  A conformer with one
  folded arm is a bad conformer however many arms it has.
"""
import inspect
import unittest

import numpy as np

from pestifer.charmmff.make_pdb_collection import (_MAX_FOLD_REDRAWS, _MIN_CHAIN_EXTENSION,
                                                   _TRANS_CC_RISE,
                                                   _sample_and_write_mc_conformers,
                                                   acyl_chain_rows, chain_extensions,
                                                   conformer_is_folded)

Z = np.array([0.0, 0.0, 1.0])


def _straight(n_chains=2, n_carbons=8, drop_fraction=1.0):
    """Coordinates for `n_chains` chains descending `drop_fraction` of their all-trans contour."""
    names, xyz = ['P'], [[0.0, 0.0, 5.0]]
    for c in range(n_chains):
        for i in range(n_carbons):
            names.append(f'C{c + 2}{i + 1}')
            xyz.append([float(c), 0.0, -i * _TRANS_CC_RISE * drop_fraction])
    return names, np.array(xyz)


class TestChainDetection(unittest.TestCase):

    def test_two_chains_are_found_with_their_ends(self):
        names, _ = _straight(n_chains=2, n_carbons=10)
        rows = acyl_chain_rows(names)
        self.assertEqual(len(rows), 2)
        for a, t, n in rows:
            self.assertEqual(n, 10)
            self.assertLess(a, t)

    def test_a_glycerol_backbone_is_not_a_chain(self):
        """NEGATIVE CONTROL, and a bug that has already been made once here.

        `C11`-`C13` is a three-carbon glycerol backbone whose "tip" sits up by the headgroup.
        An earlier fold audit parsed it as a chain and reported PMCL1 folded when it was not.
        """
        self.assertEqual(acyl_chain_rows(['C11', 'C12', 'C13', 'P', 'N']), [])

    def test_cardiolipin_arms_are_found(self):
        names = ['P1'] + [f'C{arm}{i}' for arm in 'ABCD' for i in range(1, 17)]
        self.assertEqual(len(acyl_chain_rows(names)), 4)

    def test_min_length_is_respected(self):
        names = ['P'] + [f'C2{i}' for i in range(1, 8)]      # 7 carbons, below the floor
        self.assertEqual(acyl_chain_rows(names), [])


class TestExtensionMeasure(unittest.TestCase):

    def test_all_trans_scores_one(self):
        names, X = _straight(drop_fraction=1.0)
        self.assertTrue(all(abs(r - 1.0) < 1e-9 for r in chain_extensions(X, acyl_chain_rows(names), Z)))

    def test_a_chain_lying_in_plane_scores_zero(self):
        names, X = _straight(drop_fraction=0.0)
        self.assertTrue(all(abs(r) < 1e-9 for r in chain_extensions(X, acyl_chain_rows(names), Z)))

    def test_a_chain_that_rises_scores_negative(self):
        names, X = _straight(drop_fraction=-0.5)
        self.assertTrue(all(r < 0 for r in chain_extensions(X, acyl_chain_rows(names), Z)))

    def test_each_chain_is_measured_against_its_OWN_attachment(self):
        """Not against a shared head plane.

        Chains attach at different heights; measuring from one plane would score a chain that
        hangs off a lower attachment as short when it is perfectly extended.
        """
        names, X = _straight(n_chains=2, n_carbons=8, drop_fraction=1.0)
        rows = acyl_chain_rows(names)
        X[rows[1][0]:rows[1][1] + 1, 2] -= 6.0        # drop the whole second chain bodily
        r = chain_extensions(X, rows, Z)
        self.assertAlmostEqual(r[0], r[1], places=9,
                               msg='a bodily-translated chain changed score, so the reference is '
                                   'not its own attachment')


class TestTheFoldVerdict(unittest.TestCase):

    def test_one_folded_chain_condemns_the_conformer(self):
        """The MINIMUM over chains, never the mean -- the whole shape of the POPC__Lo defect.

        Its sn-1 extended at r~0.9 while sn-2 folded at r~0.07; a mean would read ~0.5 and pass.
        """
        names, X = _straight(n_chains=2, n_carbons=8, drop_fraction=1.0)
        rows = acyl_chain_rows(names)
        a, t, n = rows[1]
        for k, row in enumerate(range(a, t + 1)):       # fold the second chain back up
            X[row][2] = 5.0 - k * 0.05
        r = chain_extensions(X, rows, Z)
        self.assertGreater(max(r), 0.9)
        self.assertLess(min(r), 0.1)
        self.assertGreater(sum(r) / len(r), _MIN_CHAIN_EXTENSION,
                           'POSITIVE CONTROL: the mean must pass, or this proves nothing')
        self.assertTrue(conformer_is_folded(X, rows, Z))

    def test_a_fully_extended_conformer_passes(self):
        names, X = _straight(drop_fraction=1.0)
        self.assertFalse(conformer_is_folded(X, acyl_chain_rows(names), Z))

    def test_a_molecule_with_no_chains_is_not_condemned(self):
        """A sterol or anything else without acyl chains must not be rejected by a measure that
        does not apply to it."""
        self.assertFalse(conformer_is_folded(np.zeros((3, 3)), [], Z))


class TestTheThresholdAgainstMeasuredValues(unittest.TestCase):
    """Real per-chain minima from the shipped collections, so the cut is pinned to evidence."""

    # Membership here is set by an INDEPENDENT geometric classifier -- a chain is folded when its
    # tip both drops less than a quarter of its contour AND stays radially close, i.e. it came
    # back up rather than splaying sideways -- not by whether this threshold happens to catch it.
    # LNACL2's 0.24 arm was mislabelled sound while writing these tests, on a hunch that
    # four-chain lipids were being punished by order statistics; the classifier says fold
    # (drop 5.09 A against a 5.40 A quarter-contour, radial 3.25 A), and the per-chain rates say
    # the order-statistics worry was unfounded.
    SHIPPED_FOLDS = {'SOPG': 0.02, 'POPC__Lo (v3.25.1 defect)': 0.07, 'DEPA': 0.07,
                     'POPI24': 0.11, 'SDPS': 0.11, 'DEPC': 0.12, 'DEPE': 0.14,
                     'DNPE': 0.16, 'DXPC': 0.21, 'LNACL2 arm': 0.24}
    SHIPPED_SOUND = {'DMPI2A': 0.26, 'POPE': 0.27, 'DCPC': 0.27,
                     'SAPI24': 0.42, 'DMPC': 0.61, 'POPC__Lo (repaired)': 0.76,
                     'POPC__Lo (pre-regression)': 0.79}

    def test_every_measured_fold_is_rejected(self):
        for name, r in self.SHIPPED_FOLDS.items():
            self.assertLess(r, _MIN_CHAIN_EXTENSION, f'{name} (r={r}) would now be accepted')

    def test_every_measured_sound_chain_is_accepted(self):
        for name, r in self.SHIPPED_SOUND.items():
            self.assertGreaterEqual(r, _MIN_CHAIN_EXTENSION, f'{name} (r={r}) would now be rejected')

    def test_the_cut_sits_in_the_gap_rather_than_on_an_edge(self):
        self.assertGreater(_MIN_CHAIN_EXTENSION, max(self.SHIPPED_FOLDS.values()))
        self.assertLessEqual(_MIN_CHAIN_EXTENSION, min(self.SHIPPED_SOUND.values()))


class TestTheSamplerAppliesIt(unittest.TestCase):
    """Aimed at the call site: a correct measure nobody consults is the original bug."""

    def setUp(self):
        self.src = inspect.getsource(_sample_and_write_mc_conformers)

    def test_conformers_are_filtered_before_being_written(self):
        self.assertIn('conformer_is_folded(', self.src)
        filt = self.src.index('conformer_is_folded(')
        write = self.src.index("open(f'{resid}-{f:0{digits}d}.pdb'")
        self.assertLess(filt, write, 'conformers are written before the fold filter runs')

    def test_a_shortfall_is_topped_up_by_redrawing(self):
        self.assertIn('_MAX_FOLD_REDRAWS', self.src)
        self.assertIn('redraw_seed', self.src)
        self.assertGreaterEqual(_MAX_FOLD_REDRAWS, 1)

    def test_a_redraw_keeps_the_tuned_bias(self):
        """A redraw at a different bias would move the chain order the tuner just resolved, so
        the ensemble would no longer be the phase its own metadata claims."""
        after = self.src[self.src.index('redraw_seed'):]
        call = after[after.index('run_mc('):]
        self.assertIn('torsion_bias=used_bias', call[:500])
        self.assertIn('seed=redraw_seed', call[:500])

    def test_it_refuses_rather_than_shipping_a_short_ensemble(self):
        self.assertIn('raise RuntimeError', self.src[self.src.index('_MAX_FOLD_REDRAWS'):])

    def test_it_is_a_tail_rejection_not_a_ranking(self):
        """What must NOT happen: keeping the most-extended N would bias every ensemble toward
        order and could thicken every bilayer built from the collection."""
        for banned in ('sorted(', '.sort(', 'nlargest', 'argsort'):
            self.assertNotIn(banned, self.src[self.src.index('chain_rows ='):],
                             f'the filter ranks conformers ({banned}) instead of rejecting a tail')


if __name__ == '__main__':
    unittest.main()
