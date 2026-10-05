# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""An ordered conformer ensemble must not buy its chain order by folding.

v3.25.1 shipped ``POPC__Lo`` with all ten conformers hairpinned: the sn-2 chain doubled back so
its tip sat 7.7 A below its own phosphate instead of 23.5.  It cost ~4.7 A of built bilayer
thickness on ex17's patchA, and nothing caught it.

Why nothing caught it is the part worth keeping.  The trans-ordering objective is
``1/2(3cos^2(theta) - 1)``, invariant under ``theta -> 180 - theta``, so an inverted chain scores
EXACTLY as well as an extended one -- folding is a free way to hit the order target.  The axial
ceiling is one plane at the lowest head reference atom, so a hairpin that comes back to just under
that plane carries zero axial excess.  And the fold audit keyed on "tip above the phosphate", which
these conformers never are.  Three separate guards, all blind to the same shape.

The control that *was* available: ``_bisect_to_order`` always evaluates the lower bound (bias 0)
first, so the same residue's unbiased ensemble is already in hand.  Ordering straightens chains, so
an ordered ensemble cannot be dramatically shorter than the fluid one it came from.

Measured 2026-10-05: of five MC seeds for POPC's Lo ensemble, only seed 0 collapsed -- and seed 0
is the one that shipped.  The other four reached the same chain order (0.294-0.308) at a trans bias
of 2.5-5.0, where seed 0 chased it to 20.0 by folding.  So this is a frozen bad draw, which is
exactly what a reseed fixes.
"""
import inspect
import unittest

from pestifer.charmmff.make_pdb_collection import (_LO_COLLAPSE_FRACTION, _MAX_ORDER_RESEEDS,
                                                   _sample_and_write_mc_conformers,
                                                   ensemble_collapsed)


class TestCollapseDetection(unittest.TestCase):

    def test_the_shipped_defect_is_caught(self):
        """The real numbers: POPC__Lo 2.56 A tuned against 15.33 A unbiased."""
        self.assertTrue(ensemble_collapsed(2.56, 15.33))

    def test_a_healthy_ordered_ensemble_is_not_flagged(self):
        """POSITIVE CONTROL: the same residue at seed 1 -- ordering makes it LONGER, not shorter."""
        self.assertFalse(ensemble_collapsed(17.24, 15.33))

    def test_ordinary_sampling_noise_is_not_flagged(self):
        """The margin this threshold has to leave.

        Across 174 base/``__Lo`` pairs in three successive releases, benign cases where Lo lands
        slightly under Ld bottom out at a ratio of 0.94 (TMCL1, 9.91 vs 10.52).  Flagging those
        would make the guard fire on almost every release and be switched off.
        """
        self.assertFalse(ensemble_collapsed(9.91, 10.52))
        self.assertFalse(ensemble_collapsed(10.33, 10.50))   # BSM, ratio 0.98

    def test_the_other_known_collapse_is_caught(self):
        """DMPG__Lo in the release before: 1.72 A against 7.21 A, ratio 0.24.

        Two real collapses at 0.17 and 0.24, worst benign case at 0.94 -- so the threshold is
        calibrated against both sides, not fitted to the one defect that prompted it.
        """
        self.assertTrue(ensemble_collapsed(1.72, 7.21))

    def test_the_threshold_sits_between_them_with_margin(self):
        self.assertGreater(_LO_COLLAPSE_FRACTION, 0.24 * 1.5,
                           'threshold too tight to catch the known collapses with margin')
        self.assertLess(_LO_COLLAPSE_FRACTION, 0.94 * 0.75,
                        'threshold close enough to benign noise to fire spuriously')

    def test_no_baseline_means_no_verdict(self):
        """NEGATIVE CONTROL: a missing control must not become a loud wrong failure.

        The caller raises when this returns True on its last attempt, so answering True on an
        unmeasurable baseline would refuse to write conformer sets for any residue whose extent
        cannot be measured -- trading a silent wrong answer for a noisy one.
        """
        self.assertFalse(ensemble_collapsed(0.0, 0.0))
        self.assertFalse(ensemble_collapsed(5.0, 0.0))
        self.assertFalse(ensemble_collapsed(float('nan'), 15.0))


class TestTheSamplerActuallyUsesIt(unittest.TestCase):
    """The predicate being right is not evidence the sampler consults it.

    This is the failure shape recorded in CLAUDE.md: a conformer-cache test that passed against
    the helper while the call site no longer called it.  These are aimed at the call site.
    """

    def setUp(self):
        self.src = inspect.getsource(_sample_and_write_mc_conformers)

    def test_the_sampler_calls_the_predicate(self):
        self.assertIn('ensemble_collapsed(', self.src,
                      'the collapse check is no longer consulted by the sampler')

    def test_the_control_is_the_unbiased_ensemble_of_the_same_residue(self):
        """Not a constant, and not another residue: the bias-0 probe of this very molecule.

        That is what makes the check free and self-calibrating -- a short lipid and a long one
        are each compared against themselves.
        """
        line = next(ln for ln in self.src.splitlines() if 'ensemble_collapsed(' in ln)
        self.assertIn('ordered_extent', line)
        self.assertIn('unbiased_extent', line)
        self.assertIn('probed[0.0]', self.src,
                      'the unbiased control is not taken from the bias-0 probe')

    def test_a_collapse_reseeds_rather_than_shipping(self):
        self.assertIn('_MAX_ORDER_RESEEDS', self.src)
        self.assertIn('used_seed += 1', self.src, 'a detected collapse does not reseed')

    def test_exhausting_the_reseeds_raises(self):
        """It must fail loudly, not fall through and write the folded set anyway."""
        self.assertIn('raise RuntimeError', self.src)
        self.assertGreaterEqual(_MAX_ORDER_RESEEDS, 1)

    def test_the_seed_actually_used_is_what_gets_recorded(self):
        """Provenance has to reproduce the entry, not the request.

        A reseeded entry recorded under the requested seed would claim to be reproducible and
        not be -- the same class of wrong as the docstring that described a reverted design.
        """
        self.assertIn("'mc_seed': int(used_seed)", self.src)
        caller = inspect.getsource(
            __import__('pestifer.charmmff.make_pdb_collection', fromlist=['do_psfgen']).do_psfgen)
        self.assertIn("mc_seed = mc_result.get('mc_seed', mc_seed)", caller,
                      'do_psfgen records the requested seed, not the one the sampler used')


if __name__ == '__main__':
    unittest.main()
