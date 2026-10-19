# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A sterol ring speared by an acyl tail is un-threaded by rotation, not by deleting a molecule.

``ring_check``'s rotation engine was reachable only by aromatic side chains, glycan-speared rigid
rings, and protein-speared glycan rings. A lipid ring speared by a lipid bond matched none of
them, so **every membrane piercing went straight to deletion without one rotation attempted** —
an ex17 build deleted 7-11 lipids per run that way. At one or two piercings that is tolerable; at
ten it changes the composition of the membrane being built, and the regenerated conformer
collection roughly doubled the count (quilt piercings 5/3/7 under 3.26.0, 10/13/8 under 3.27.0).

Measured on a real pre-``ring_check`` ex17 quilt, all 15 piercings were the same motif — a
cholesterol ring speared by an acyl tail — and **15 of 15 cleared by bending the tail about an
upstream chain dihedral. Zero deletions were needed.**

Two moves are wired, cheapest-first by how much they disturb:

1. **Tail dihedral** — delegates to ``try_glycan_pendant``, which is chemistry-agnostic by its own
   docstring: it rigidly rotates a graph branch. The branch is a few tail carbons.
2. **Sterol long-axis spin** — turns the speared molecule in place. It cleared **0 of 15** here,
   because a fused ring plate is roughly coaxial with its long axis, so spinning moves the ring
   *atoms* much more than it moves the *hole*. It is kept as the remedy for a piercing bond with
   no usable upstream hinge, and is so far unexercised on real data. Recorded as such rather than
   presented as a working second mechanism.

Deletion remains, as the last resort, and now says rotation was tried.
"""
import inspect
import unittest

import numpy as np

from pestifer.psfutil.ring_resolve import _SPIN_DEGREES, try_sterol_spin
from pestifer.tasks.ringcheck import RingCheckTask


def _p(piercee_type, piercer_type, **kw):
    pe = dict(segtype=piercee_type, resname=kw.get('pe_resname', 'CHL1'),
              segname='L1', resid='1', ring_serials=kw.get('ring', [1, 2, 3, 4, 5, 6]))
    pr = dict(segtype=piercer_type, resname=kw.get('pr_resname', 'SOPE'),
              segname='L1', resid='2', bond_serials=kw.get('bond', [10, 11]))
    return {'piercee': pe, 'piercer': pr}


class TestLipidPiercingsReachTheRotationEngine(unittest.TestCase):

    def test_a_sterol_speared_by_an_acyl_tail_is_rotatable(self):
        """The whole defect in one assertion: this used to be False."""
        self.assertTrue(RingCheckTask._rotatable(_p('lipid', 'lipid')))

    def test_it_still_needs_something_to_rotate(self):
        p = _p('lipid', 'lipid')
        p['piercer']['bond_serials'] = []
        p['piercee']['ring_serials'] = []
        self.assertFalse(RingCheckTask._rotatable(p))

    def test_the_other_kinds_are_unaffected(self):
        """NEGATIVE CONTROL: widening the predicate must not disturb the three it already had."""
        self.assertTrue(RingCheckTask._rotatable(
            dict(piercee=dict(segtype='protein', resname='PHE'), piercer=dict(segtype='lipid'))))
        self.assertFalse(RingCheckTask._rotatable(
            dict(piercee=dict(segtype='protein', resname='ALA'),
                 piercer=dict(segtype='lipid', bond_serials=[1, 2]))))

    def test_a_lipid_piercing_is_NOT_a_water_or_ion_piercing(self):
        for other in ('water', 'ion', 'unknown'):
            self.assertFalse(RingCheckTask._rotatable(_p(other, 'lipid')))
            self.assertFalse(RingCheckTask._rotatable(_p('lipid', other)))


class TestTheDispatchOrder(unittest.TestCase):
    """Tail first, spin second. The tail move disturbs a few carbons; the spin turns a whole
    molecule, so trying it first would pay more collateral motion for the same result."""

    def setUp(self):
        self.src = inspect.getsource(RingCheckTask._resolve_by_rotation)

    def test_both_lipid_moves_are_wired(self):
        self.assertIn("mover_kind='lipid tail'", self.src)
        self.assertIn('try_sterol_spin', self.src)

    def test_the_tail_move_is_tried_first(self):
        self.assertLess(self.src.index("mover_kind='lipid tail'"),
                        self.src.index('try_sterol_spin'))

    def test_the_spin_only_runs_if_the_tail_move_failed(self):
        seg = self.src[self.src.index("mover_kind='lipid tail'"):self.src.index('try_sterol_spin')]
        self.assertIn('if cand is None', seg)


class TestTheResolutionLoop(unittest.TestCase):
    """Clearing a piercing can create another, and a molecule can be threaded twice -- which one
    rigid rotation about one hinge cannot undo.  So resolution iterates.

    These replace three source-text assertions that pinned the single-pass shape: that the lipid
    branch did not raise, that a non-lipid one did, and that the copy was guarded.  The loop
    moved all three, and the behaviour they were protecting is covered directly by
    ``test_ringcheck_task.TestAnUnresolvableLipidFallsThroughInsteadOfCrashing`` and
    ``TestFatalPiercings``.  What is left here is what only the source can show: that there IS a
    bounded loop and that nothing is copied when nothing was written.
    """

    def setUp(self):
        self.src = inspect.getsource(RingCheckTask.do)

    def test_resolution_is_a_bounded_loop(self):
        self.assertIn('_MAX_RESOLVE_CYCLES', self.src, 'resolution no longer iterates')
        from pestifer.tasks.ringcheck import _MAX_RESOLVE_CYCLES
        self.assertGreater(_MAX_RESOLVE_CYCLES, 1, 'one pass is not a loop')
        self.assertLess(_MAX_RESOLVE_CYCLES, 20, 'an unbounded retry on a structure that cannot '
                                                 'be fixed is a hang, not a resolution')

    def test_nothing_is_copied_when_no_pose_was_written(self):
        """`_resolve_by_rotation` returns None for the pose when it clears nothing, and
        `shutil.copy(None, ...)` raised `TypeError: stat: path should be string ... not
        NoneType` on the first build that reached it."""
        i = self.src.index('shutil.copy(fixed_pdb')
        before = self.src[:i]
        self.assertIn('if fixed_pdb is None', before,
                      'the copy is reachable with fixed_pdb None again')

    def test_the_loop_stops_when_a_pass_makes_no_progress(self):
        """Otherwise a structure whose piercings cannot be cleared burns every cycle doing the
        same work."""
        self.assertIn('no net', self.src, 'no no-progress exit')


class TestTheSterolSpin(unittest.TestCase):

    def test_it_rotates_about_the_long_axis_through_the_centroid(self):
        """Position and membrane-normal alignment must be preserved: a sterol displaced
        laterally would trade a piercing for a clash."""
        src = inspect.getsource(try_sterol_spin)
        self.assertIn('centroid', src)
        self.assertIn('svd', src)

    def test_it_declines_a_molecule_too_small_to_define_an_axis(self):
        class _Ck:
            _atoms = []
            _row_of_serial = {}
        self.assertIsNone(try_sterol_spin(_Ck(), 'x.pdb', np.zeros((3, 3)), None,
                                          dict(segname='L1', resid='1', resname='X',
                                               segtype='lipid'), ('L1', '1'), 'o.pdb'))

    def test_the_sweep_leads_with_the_largest_turns(self):
        """The acceptance test is binary -- does it un-thread -- so a fine sweep only costs time."""
        self.assertEqual(_SPIN_DEGREES[0], 180)
        self.assertEqual(len(set(_SPIN_DEGREES)), len(_SPIN_DEGREES))
        self.assertTrue(all(0 < d < 360 for d in _SPIN_DEGREES))


if __name__ == '__main__':
    unittest.main()


class TestPendantAxesIsScopedToTheMOLECULE(unittest.TestCase):
    """``pendant_axes`` built its graph from the whole segname and then copied that graph once
    per bridge. A glycan segment holds one molecule, so that was fine and never noticed.

    A membrane puts every lipid in one segment. Measured on an ex17 quilt: segment ``L1`` holds
    **56,535 atoms in 538 molecules with 52,173 bridges**, so the loop would copy a 56,535-node
    graph 52,173 times — about three hours per piercing, ten to thirteen times per build. It
    would have presented as a hang, not as a bug.

    Scoped to the molecule carrying the piercing bond it takes milliseconds, and the answer is
    identical: a hinge in a different molecule cannot move this bond. This test pins the
    correctness half of that, which is the durable half — a timing test would be flaky and would
    not say *why* the scope is right.
    """

    def test_hinges_come_only_from_the_piercers_own_molecule(self):
        import networkx as nx
        from unittest import mock
        from pestifer.psfutil.psfring import RingChecker

        # two disjoint 6-atom chains in ONE segment, as lipids are
        G = nx.Graph()
        G.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 5), (5, 6)])       # molecule A
        G.add_edges_from([(11, 12), (12, 13), (13, 14), (14, 15), (15, 16)])  # molecule B
        rc = mock.Mock(spec=RingChecker)
        rc._atoms = [mock.Mock(serial=s, segname='L1') for s in list(G.nodes())]
        rc.topol = mock.Mock(G=G)
        rc._row_of_serial = {s: i for i, s in enumerate(sorted(G.nodes()))}

        axes = RingChecker.pendant_axes(rc, [3, 4], 'L1')
        self.assertTrue(axes, 'no hinge found for a bond that has upstream bridges')
        molecule_a = set(G.nodes()) - {11, 12, 13, 14, 15, 16}
        for ax in axes:
            self.assertIn(ax['i_ser'], molecule_a,
                          'a hinge was taken from a different molecule in the same segment')
            self.assertIn(ax['j_ser'], molecule_a)
            moved = {s for s, r in rc._row_of_serial.items() if r in set(ax['prows'].tolist())}
            self.assertTrue(moved <= molecule_a,
                            f'the moving branch escapes the piercer\'s molecule: {moved}')


class TestASterolSpearedByAProteinSideChain(unittest.TestCase):
    """The one combination with no entry in ``_rotatable``: a LIPID ring pierced by a PROTEIN
    bond.

    It fell through every branch to the lipid deletion path, so ``delete: piercee`` removed the
    **cholesterol** and left the side chain doing the threading exactly where it was -- the
    sterol depletion ``d89cc518`` eliminated, re-entering through a different door, and not even
    a fix for the piercing.

    Nothing had ever reached this branch because no example puts ``protein`` in ``ring_check``'s
    ``segtypes``.  It was found on the replica of the 2026-10-16 example-17 arm that died
    post-embed: ``CHL1 159 ring <- LEU 702 bond``, with a second piercing ``PHE 699 ring <-
    SOPE 119 bond`` on the same structure.  The dying atom was ``LEU 702 CB``.
    """

    def test_a_sterol_speared_by_a_protein_bond_is_rotatable(self):
        """The defect in one assertion: this used to be False, and False meant delete a sterol."""
        self.assertTrue(RingCheckTask._rotatable(_p('lipid', 'protein', pr_resname='LEU')))

    def test_it_still_needs_something_to_rotate(self):
        p = _p('lipid', 'protein', pr_resname='LEU')
        p['piercer']['bond_serials'] = []
        p['piercee']['ring_serials'] = []
        self.assertFalse(RingCheckTask._rotatable(p))

    def test_the_reverse_pairing_was_already_handled_and_still_is(self):
        """An AROMATIC protein ring speared by a lipid bond is case 1 and must be untouched."""
        self.assertTrue(RingCheckTask._rotatable(
            dict(piercee=dict(segtype='protein', resname='PHE'),
                 piercer=dict(segtype='lipid', bond_serials=[1, 2]))))

    def test_a_non_aromatic_protein_ring_speared_by_a_lipid_is_still_not_rotatable(self):
        """NEGATIVE CONTROL. Widening for lipid<-protein must not also widen protein<-lipid:
        there is no ring to spin on a LEU, so it has no remedy and must not claim one."""
        self.assertFalse(RingCheckTask._rotatable(
            dict(piercee=dict(segtype='protein', resname='LEU'),
                 piercer=dict(segtype='lipid', bond_serials=[1, 2]))))

    def test_water_and_ion_piercers_are_still_rejected(self):
        for other in ('water', 'ion', 'unknown'):
            self.assertFalse(RingCheckTask._rotatable(_p('lipid', other)))


class TestTheProteinPiercerDispatchOrder(unittest.TestCase):
    """Spin the sterol first, disturb the protein second -- case 3's rule applied here: move the
    model-built partner before the deposited one, and in a membrane build the lipid is the
    grid-placed one.

    ``try_sterol_spin`` is honest about its own limit (it clears a bond near the ring rim and
    returns ``None`` for one through the centre), so the fall-through costs one attempt and the
    side-chain rotation is the move that actually pulls the bond out of the ring.
    """

    def setUp(self):
        self.src = inspect.getsource(RingCheckTask._resolve_by_rotation)

    def test_the_protein_piercer_branch_exists(self):
        self.assertIn("piercer['segtype'] == 'protein'", self.src)
        self.assertIn('_try_piercer_sidechain', self.src)

    def test_the_sterol_spin_is_attempted_before_the_side_chain(self):
        branch = self.src[self.src.index("# 5. lipid ring speared by a PROTEIN bond"):]
        self.assertLess(branch.index('try_sterol_spin'), branch.index('_try_piercer_sidechain'),
                        'the side chain is being disturbed before the sterol is spun')

    def test_the_branch_is_guarded_by_cand_is_none(self):
        """Every earlier remedy must get its turn first; an unguarded branch would re-resolve a
        piercing another case had already cleared."""
        branch = self.src[self.src.index("# 5. lipid ring speared by a PROTEIN bond"):]
        guard = branch[branch.index('if ('):branch.index('cand = try_sterol_spin')]
        self.assertIn('cand is None', guard)
