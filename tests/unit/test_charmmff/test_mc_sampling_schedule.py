# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""``run_mc`` must collect a conformer on the step schedule, not only when a move is accepted.

The sample-collection test used to sit at the BOTTOM of the proposal loop, after
``n_accept += 1`` -- and every rejection path reaches that point by ``continue``, so a checkpoint
fired only if the proposal that happened to land on it was accepted.  At the acceptance rates a
confined, overlap-checked lipid actually runs at, most checkpoints never fired, and the shortfall
was padded with copies of the final conformer by a line whose comment called it "the unlikely
event decorrelation math left us short".

It was not unlikely; it was the normal path.  The shipped 3.26.0 collection held **1416 distinct
conformers out of 2690 (53%)**, no entry reached ten distinct, and ``DMPI33__Lo`` was a single
conformer repeated ten times.  POPC shipped four distinct of ten, one draw repeated seven times
-- which is what a rendered gallery page made obvious at a glance.

Rejecting a move does not mean there is nothing to record: in Metropolis MC the chain stays in
its current state and that state is re-counted.  Sampling belongs on the step schedule.

Regenerating POPC with the fix, same seed and settings: 10 distinct of 10, and the four
previously-shipped conformers are all still among them -- the fix adds the six that padding was
swallowing rather than changing where the sampler goes.
"""
import inspect
import unittest

import networkx as nx
import numpy as np

from pestifer.charmmff.athermal_mc import (MoleculeMC, RotatableBond, build_exclusions,
                                           moving_set, run_mc)


def _trans_chain(n=16, bond=1.53, angle_deg=112.0):
    half = np.deg2rad(angle_deg) / 2.0
    dz, dx = bond * np.sin(half), bond * np.cos(half)
    coords = np.zeros((n, 3))
    for i in range(1, n):
        coords[i, 2] = coords[i - 1, 2] + dz
        coords[i, 0] = coords[i - 1, 0] + (dx if i % 2 else -dx)
    return coords, nx.path_graph(n), np.full(n, 1.7 * 2.0 ** (-1.0 / 6.0))


def _chain_mol(n=16, cylinder_radius=float('inf')):
    coords, g, radii = _trans_chain(n)
    rot = [RotatableBond(a=i, b=i + 1, moving=moving_set(g, i, i + 1)) for i in range(1, n - 1)]
    return MoleculeMC(coords=coords, radii=radii, rotatable=rot,
                      exclusions=build_exclusions(g, order=2),
                      axis_point=coords[0], axis_dir=np.array([0.0, 0.0, 1.0]),
                      cylinder_radius=cylinder_radius)


def _distinct(samples):
    return len({np.round(s, 6).tobytes() for s in samples})


class TestEveryCheckpointIsCollected(unittest.TestCase):

    def test_a_free_chain_returns_all_distinct_conformers(self):
        s = run_mc(_chain_mol(), nsamples=10, n_equil=400, n_decorr=120, seed=5)
        self.assertEqual(len(s), 10)
        self.assertEqual(_distinct(s), 10)

    def test_a_tightly_confined_chain_still_returns_all_distinct(self):
        """The case that was failing in production.

        Confinement drives the acceptance rate down, which is exactly when the old code lost
        checkpoints: the tighter the cylinder, the more duplicates it produced.  Diversity must
        not depend on how often a proposal is accepted.
        """
        s = run_mc(_chain_mol(cylinder_radius=3.0), nsamples=10, n_equil=400, n_decorr=120, seed=5)
        self.assertEqual(len(s), 10)
        self.assertEqual(_distinct(s), 10)

    def test_no_padding_warning_on_an_ordinary_run(self):
        """The padding path must be unreachable in normal operation, not merely quiet."""
        with self.assertNoLogs('pestifer.charmmff.athermal_mc', level='WARNING'):
            run_mc(_chain_mol(cylinder_radius=3.0), nsamples=8, n_equil=300, n_decorr=100, seed=2)

    def test_sample_count_is_independent_of_the_acceptance_rate(self):
        """POSITIVE CONTROL that the two arms really do differ in acceptance.

        If a tight cylinder did not actually lower acceptance, the test above would be comparing
        two identical situations and would prove nothing.
        """
        loose, tight = _chain_mol(), _chain_mol(cylinder_radius=2.5)
        for mol in (loose, tight):
            s = run_mc(mol, nsamples=6, n_equil=300, n_decorr=100, seed=9)
            self.assertEqual(len(s), 6)
        # a chain that cannot move at all is the degenerate end: identical samples are then
        # CORRECT, because the state really did not change
        frozen = _chain_mol(cylinder_radius=0.01)
        s = run_mc(frozen, nsamples=5, n_equil=200, n_decorr=50, seed=1)
        self.assertEqual(len(s), 5)


class TestTheScheduleIsNotGatedOnAcceptance(unittest.TestCase):
    """Aimed at the mechanism, because the symptom only shows at production step counts."""

    @staticmethod
    def _code_only(text):
        """Source with comments and docstring text stripped.

        Structural assertions must read the CODE.  A first version of these tests matched its own
        explanatory comment -- which quotes `n_accept += 1` and `samples.append` verbatim -- and
        located the wrong line.  That is the same shape as a verification script that passed
        because a comment contained the string it was grepping for.
        """
        out = []
        for line in text.splitlines():
            stripped = line.split('#', 1)[0]
            if stripped.strip():
                out.append(stripped)
        return '\n'.join(out)

    def setUp(self):
        self.src = self._code_only(inspect.getsource(run_mc))
        self.body = self.src[self.src.index('for step in range('):]

    def test_collection_does_not_follow_the_accept_counter(self):
        """The literal defect: `samples.append` after `n_accept += 1` and still inside the loop.

        Bounded by indentation rather than by a blank line: with comments stripped there are no
        blank lines left, so a naive split runs past the loop and swallows the post-loop
        checkpoint and the padding, both of which legitimately append.
        """
        lines = self.body.splitlines()
        start = next(i for i, ln in enumerate(lines) if 'n_accept += 1' in ln)
        tail = []
        for ln in lines[start + 1:]:
            if len(ln) - len(ln.lstrip()) <= 4:     # dedented out of the for-loop body
                break
            tail.append(ln)
        self.assertNotIn('samples.append', '\n'.join(tail),
                         'sampling is still gated on a move being accepted')

    def test_the_checkpoint_precedes_the_proposal(self):
        head = self.body[:self.body.index('n_accept += 1')]
        self.assertIn('samples.append', head,
                      'no checkpoint runs before the proposal logic, so rejections skip it')

    def test_a_final_checkpoint_runs_after_the_loop(self):
        """Taking checkpoints at the top of the NEXT iteration leaves the last one unreachable
        from inside the loop; without this the ensemble is one conformer short every time."""
        tail = self.src[self.src.index('n_accept += 1'):]
        self.assertIn('total_proposals - n_equil', tail)

    def test_padding_survives_but_announces_itself(self):
        """It stays as a safety net -- but silently absorbing missed checkpoints is how the
        duplicates reached the collection, so if it fires again that must be visible."""
        self.assertIn('while len(samples) < nsamples', self.src)
        guard = self.src[self.src.index('if len(samples) < nsamples'):]
        self.assertIn('logger.warning', guard[:guard.index('while len(samples)')])


class TestDeterminismIsPreserved(unittest.TestCase):
    """What must NOT change: the sampler's trajectory is still seed-reproducible."""

    def test_same_seed_same_ensemble(self):
        a = run_mc(_chain_mol(), nsamples=6, n_equil=300, n_decorr=100, seed=4)
        b = run_mc(_chain_mol(), nsamples=6, n_equil=300, n_decorr=100, seed=4)
        for x, y in zip(a, b):
            np.testing.assert_allclose(x, y)

    def test_different_seeds_differ(self):
        a = run_mc(_chain_mol(), nsamples=6, n_equil=300, n_decorr=100, seed=4)
        b = run_mc(_chain_mol(), nsamples=6, n_equil=300, n_decorr=100, seed=5)
        self.assertFalse(np.allclose(a[-1], b[-1]))


if __name__ == '__main__':
    unittest.main()
