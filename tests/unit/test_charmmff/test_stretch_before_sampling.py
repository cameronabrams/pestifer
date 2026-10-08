# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A multi-arm lipid is stretched before the MC sampler melts it; a two-chain lipid is not.

The colvar-driven stretch pushes tails into a rod, which is exactly what the athermal-MC sampler
exists to avoid, so it was skipped for ``sampler='mc'`` outright.  Measured 2026-10-08, skipping
it is the worse trade on a crowded molecule.

A cardiolipin's four chains are built from internal coordinates into a coil around one backbone,
and torsional MC with overlap rejection opens three of them and leaves the fourth closed.  Across
the six CL2-type cardiolipins the same arm sits in a band of roughly **0.24-0.28** of its
all-trans contour while its three siblings span 0.4-0.9 -- about ten times narrower, and
straddling the 0.25 fold threshold, so which side an entry falls on is near-arbitrary.  LNACL2
lands below and could not be generated at all; TLCL2 clears it by a hundredth.

Running the stretch first lifts that arm to 0.77 (LNACL2: A 0.82, B 0.77, C 0.67, D 0.62), so
nothing about the molecule prevents it -- only where the sampler starts.  Melting a rod back into
a fluid ensemble is the one thing this sampler is known to do well.

Two-chain lipids are left alone on purpose: their worst arm already has a median of 0.71 straight
from the IC build, so a rod start would buy nothing and risk the over-ordering the MC path exists
to prevent.

Recorded here because the diagnosis was wrong twice before it was right. The arm was first called
"threaded" -- it is not; rendering it with its siblings hidden shows a plain coil passing through
nothing. And TLCL2 was offered as a healthy control when its own arm sits at 0.25: the comparator
was the same phenomenon on the other side of the line.
"""
import inspect
import unittest

import numpy as np

from pestifer.charmmff.make_pdb_collection import (_MAX_CHAINS_WITHOUT_STRETCH,
                                                   _MIN_IC_END_TO_END, _TRANS_CC_RISE, do_psfgen,
                                                   min_end_to_end_over_contour,
                                                   stretch_before_sampling)

# IC-build end-to-end / contour, measured on the shipped structures.  These are what the cut is
# calibrated against, and the separation is the whole argument for the rule.
CL2 = 0.28          # every CL2-type cardiolipin, and LNACL2 -- all fixed by stretching
PMCL2_EE = 0.36     # stopped generating entirely when it was stretched
CL1 = 0.44          # already healthy; stretching regressed several
COILED_TWO = 0.28   # DDOPC, SDPC, PDOPE and nine others -- coiled, and opened fine unaided


class TestWhoGetsStretched(unittest.TestCase):

    def test_the_md_sampler_always_stretches(self):
        """Unchanged behaviour: its ensemble comes out of the stretched state."""
        for n in (1, 2, 3, 4, 6):
            self.assertTrue(stretch_before_sampling('md', n, 0.9), f'{n} chains')

    def test_a_crowded_lipid_built_coiled_is_stretched(self):
        self.assertTrue(stretch_before_sampling('mc', 4, CL2))

    def test_a_crowded_lipid_built_open_is_NOT_stretched(self):
        """Crowding alone was the first rule, and it was wrong.

        Stretching every cardiolipin fixed the CL2 family but regressed the CL1 entries and
        stopped PMCL1/PMCL2 generating at all: 12 improved against 8 worse or failing.
        """
        self.assertFalse(stretch_before_sampling('mc', 4, CL1))
        self.assertFalse(stretch_before_sampling('mc', 4, PMCL2_EE))

    def test_a_coiled_TWO_chain_lipid_is_NOT_stretched(self):
        """A coiled start alone was the other candidate rule, and it over-selects.

        Twelve two-chain entries are built just as coiled as the cardiolipins -- DDOPC, SDPC,
        PDOPE and the rest -- and the sampler opens every one of them unaided (worst-arm medians
        0.48 to 0.83). Stretching them would be churn on working entries.
        """
        self.assertFalse(stretch_before_sampling('mc', 2, COILED_TWO))

    def test_three_chains_already_counts_as_crowded(self):
        """The boundary is 'more than two', not 'four'.

        Nothing in the collection has exactly three acyl chains today, but the reason is crowding
        around one backbone, not the number four, so a triacylglycerol added later is covered
        without anyone revisiting this.
        """
        self.assertTrue(stretch_before_sampling('mc', 3, CL2))

    def test_the_cut_separates_the_measured_groups_with_margin(self):
        """Calibrated against both sides, not fitted to the entries that motivated it."""
        self.assertGreater(_MIN_IC_END_TO_END, CL2)
        self.assertLess(_MIN_IC_END_TO_END, PMCL2_EE)
        self.assertLess(_MIN_IC_END_TO_END, CL1)
        self.assertEqual(_MAX_CHAINS_WITHOUT_STRETCH, 2)


class TestTheEndToEndMeasure(unittest.TestCase):

    @staticmethod
    def _chain(n=10, spread=1.0):
        """One chain of n carbons laid along x at `spread` of the all-trans rise per bond."""
        names = ['P'] + [f'C2{i + 1}' for i in range(n)]
        xyz = [[0.0, 0.0, 5.0]] + [[i * _TRANS_CC_RISE * spread, 0.0, 0.0] for i in range(n)]
        return names, np.array(xyz)

    def test_a_fully_extended_chain_scores_one(self):
        from pestifer.charmmff.make_pdb_collection import acyl_chain_rows
        names, xyz = self._chain(spread=1.0)
        self.assertAlmostEqual(min_end_to_end_over_contour(xyz, acyl_chain_rows(names)), 1.0, 9)

    def test_a_coiled_chain_scores_low(self):
        from pestifer.charmmff.make_pdb_collection import acyl_chain_rows
        names, xyz = self._chain(spread=0.2)
        self.assertLess(min_end_to_end_over_contour(xyz, acyl_chain_rows(names)), _MIN_IC_END_TO_END)

    def test_it_does_not_depend_on_orientation(self):
        """NEGATIVE CONTROL, and the reason it is end-to-end rather than an axial drop.

        The measurement is taken before anything orients the molecule along the membrane normal,
        so a measure that read the z axis would be reading an arbitrary direction.
        """
        from pestifer.charmmff.make_pdb_collection import acyl_chain_rows
        names, xyz = self._chain(spread=0.8)
        rows = acyl_chain_rows(names)
        before = min_end_to_end_over_contour(xyz, rows)
        th = np.deg2rad(37.0)
        R = np.array([[np.cos(th), -np.sin(th), 0], [np.sin(th), np.cos(th), 0], [0, 0, 1]])
        self.assertAlmostEqual(min_end_to_end_over_contour(xyz @ R.T, rows), before, 9)

    def test_no_chains_is_not_treated_as_coiled(self):
        self.assertEqual(min_end_to_end_over_contour(np.zeros((3, 3)), []), float('inf'))
        self.assertFalse(stretch_before_sampling('mc', 4, float('inf')))


class TestTheCallSiteUsesIt(unittest.TestCase):
    """A correct predicate nobody consults is the original bug, twice over in this file's history."""

    @staticmethod
    def _code(text):
        return '\n'.join(ln.split('#', 1)[0] for ln in text.splitlines()
                         if ln.split('#', 1)[0].strip())

    def setUp(self):
        self.src = self._code(inspect.getsource(do_psfgen))

    def test_the_tasklist_consults_the_predicate(self):
        self.assertIn('stretch_before_sampling(', self.src,
                      'do_psfgen no longer asks whether to stretch; the skip is hard-coded again')

    def test_it_is_asked_with_both_the_chain_count_and_the_measurement(self):
        line = next(ln for ln in self.src.splitlines() if 'stretch_before_sampling(' in ln)
        self.assertIn('chain_info', line, f'no chain count in the decision: {line.strip()!r}')
        self.assertIn('ic_min_ee', line, f'no coiled-ness measurement: {line.strip()!r}')

    def test_the_measurement_is_taken_from_the_IC_BUILD(self):
        """Not from the minimized or oriented structure.

        A 500-step downhill minimize cannot open a closed chain, so the IC build is the state the
        decision is about -- and it is the only one available before the tasklist runs.
        """
        self.assertIn("-init.pdb", self.src)
        self.assertIn('min_end_to_end_over_contour(', self.src)

    def test_the_old_unconditional_skip_is_gone(self):
        """NEGATIVE CONTROL for a partial revert.

        If `sampler != 'mc'` reappears as the whole condition, the predicate can still exist and
        be exported and tested while every cardiolipin goes back to an unopened arm.
        """
        self.assertNotIn("if sampler != 'mc':", self.src)


if __name__ == '__main__':
    unittest.main()
