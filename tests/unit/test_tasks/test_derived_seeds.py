# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Lipid placement and loop closure take their seed from the build's ``--seed``.

Until 3.29.6 both had fixed seeds that ``--seed`` never reached: ``--seed`` sets ``namd.seed``
and nothing else.  Every replica of a membrane example started from the same lipid placement --
on example 17, the same 31 threadings in the lower-leaflet calibration patch, and the same doubly
threaded cholesterol stopping all three replicas.  Replicas that share a starting structure are
correlated, and a placement defect repeats in every one of them instead of showing at its real
rate.

These tests drive the real call paths -- the seed that reaches ``Bilayer.write_grid_pdb`` and
``close_one_loop`` -- not just the helper, because a correct helper that a call site does not use
is the failure this repository keeps a list of.
"""
import pathlib
import unittest
from unittest import mock

import numpy as np
import yaml

from pestifer.tasks.basetask import BaseTask
from pestifer.tasks.ligate import LigateTask
from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask
from pestifer.tasks.psfgen import PsfgenTask


def _provisions(seed):
    return {'namd_global_config': {} if seed is None else {'seed': seed}}


class TestDerivedSeed(unittest.TestCase):

    def _task(self, namd_seed):
        t = LigateTask.__new__(LigateTask)      # BaseTask is abstract; the method is inherited
        t.provisions = _provisions(namd_seed)
        return t

    def test_an_explicit_seed_wins(self):
        self.assertEqual(self._task(5).derived_seed(99), 99)

    def test_otherwise_the_build_seed_is_used(self):
        self.assertEqual(self._task(5).derived_seed(None), 5)

    def test_no_build_seed_falls_back_to_the_legacy_constant(self):
        self.assertEqual(self._task(None).derived_seed(None), BaseTask._LEGACY_SEED)

    def test_the_legacy_constant_is_namd_seeds_default(self):
        """The promise that a build run without --seed is unchanged rests on this equality: the
        old fixed seeds were 27021972, and so is namd.seed's schema default."""
        import pestifer
        spec = yaml.safe_load((pathlib.Path(pestifer.__file__).parent / 'schema' / 'base.yaml').read_text())
        namd = next(a for a in spec['attributes'] if a['name'] == 'namd')
        seed = next(a for a in namd['attributes'] if a['name'] == 'seed')
        self.assertEqual(seed['default'], BaseTask._LEGACY_SEED)


class TestPlacementSeed(unittest.TestCase):
    """Every placement goes through `grid_patch`, which decides the seed itself."""

    def _place(self, namd_seed, bilayer_seed=None):
        t = MakeMembraneSystemTask.__new__(MakeMembraneSystemTask)
        t.provisions = _provisions(namd_seed)
        t.bilayer_specs = {} if bilayer_seed is None else {'seed': bilayer_seed}
        t.outcome = {}
        t.basename = 'b'
        t.next_basename = mock.Mock()
        t.register = mock.Mock()
        patch = mock.Mock()
        patch.register_species_pdbs = []
        t.grid_patch(patch, patch_name='patchB')
        return patch.write_grid_pdb.call_args.kwargs['seed'], t.outcome

    def test_replicas_with_different_seeds_get_different_placements(self):
        a, _ = self._place(27021972)
        b, _ = self._place(27021973)
        self.assertNotEqual(a, b)
        self.assertEqual(b, 27021973)

    def test_an_explicit_bilayer_seed_pins_the_placement(self):
        a, _ = self._place(27021972, bilayer_seed=7)
        b, _ = self._place(27021973, bilayer_seed=7)
        self.assertEqual(a, b)
        self.assertEqual(a, 7)

    def test_the_seed_is_recorded(self):
        seed, outcome = self._place(42)
        self.assertEqual(outcome.get('placement_seed_patchB'), 42)


class TestLoopClosureSeed(unittest.TestCase):

    class _Stop(Exception):
        pass

    def _first_seed(self, namd_seed, ccd_seed=None):
        t = LigateTask.__new__(LigateTask)
        t.provisions = _provisions(namd_seed)
        t.outcome = {}
        t.basename = 'b'
        t.next_basename = mock.Mock()
        state = mock.Mock()
        state.pdb.name = 'in.pdb'
        t.get_current_artifact = lambda k: state
        t.base_molecule = mock.Mock()
        t.base_molecule.protein_loop_gaps.return_value = [{'segname': 'A', 'loop_resids': [10, 11],
                                                         'n_anchor_resid': 9, 'c_anchor_resid': 12}]
        seen = {}

        def close_one_loop(**kw):
            seen['seed'] = kw['seed']
            raise self._Stop()

        specs = {} if ccd_seed is None else {'seed': ccd_seed}
        with mock.patch('pestifer.psfutil.loop_ccd.loop_atoms_from_pdb', return_value=([], [], [])), \
             mock.patch('pestifer.psfutil.loop_ccd.close_one_loop', side_effect=close_one_loop):
            with self.assertRaises(self._Stop):
                t.close_loops_ccd(specs)
        return seen['seed'], t.outcome

    def test_the_closure_takes_the_build_seed(self):
        seed, outcome = self._first_seed(27021974)
        self.assertEqual(seed, 27021974, 'the first loop must close from the build seed + 0')
        self.assertEqual(outcome.get('ccd_seed'), 27021974)

    def test_an_explicit_ccd_seed_pins_the_closure(self):
        self.assertEqual(self._first_seed(27021974, ccd_seed=3)[0], 3)


class TestDeclashSeed(unittest.TestCase):
    """Tail modeling and the loop and pendant declashers draw from `PsfgenTask.declash_seed`.

    Each test drives a real call site and captures what reaches the modeler or the declasher, so
    a site that went back to the class constant would fail here even with `declash_seed` intact.
    """

    class _Stop(Exception):
        pass

    def _task(self, namd_seed):
        t = PsfgenTask.__new__(PsfgenTask)
        t.provisions = _provisions(namd_seed)
        t.outcome = {}
        t.basename = 'b'
        t.next_basename = mock.Mock()
        return t

    class _CM:
        def __init__(self, *a):
            self.coords = np.zeros((2, 3))
            self._mass = np.array([12.0, 12.0])
            self._names = np.array(['CA', 'CB'])
            self.atoms = mock.Mock(data=[mock.Mock(segname='A', resid='1')] * 2)

        def _residue_index(self):
            return np.zeros(2, dtype=int)

        def residue_of_segname(self, *a):
            return 0

        def write_pdb(self, *a):
            pass

    def _first_draw(self, rng):
        return int(rng.integers(1 << 30))

    def _pendant_rng(self, namd_seed):
        seen = {}

        def declash_pendant(*args):
            seen['rng'] = args[-1]
        state = mock.Mock()
        with mock.patch('pestifer.molecule.coordmanip.CoordManipulator', self._CM), \
             mock.patch('pestifer.psfutil.declash.declash_pendant', side_effect=declash_pendant):
            self._task(namd_seed)._run_pendant_declash(state, [(np.array([0]), [], [])], 1, 2.0)
        return self._first_draw(seen['rng'])

    def _loop_rng(self, namd_seed):
        seen = {}

        def declash_loop(*args):
            seen['rng'] = args[-1]
        state = mock.Mock()
        with mock.patch('pestifer.molecule.coordmanip.CoordManipulator', self._CM), \
             mock.patch('pestifer.psfutil.declash.declash_loop', side_effect=declash_loop):
            self._task(namd_seed)._run_loop_declash(state, [('A', ['1'])], 1)
        return self._first_draw(seen['rng'])

    def _tail_seed(self, namd_seed):
        t = self._task(namd_seed)
        t.specs = {'source': {'sequence': {'loops': {'declash': {}}}}}
        t.base_molecule = mock.Mock()
        t.base_molecule.protein_terminal_tails.return_value = [
            {'segname': 'A', 'tail_resids': [1, 2], 'anchor_resid': 3, 'end': 'N'}]
        t.base_molecule.protein_loop_gaps.return_value = []
        t._crotation_targeted_resids = lambda: {}
        state = mock.Mock()
        state.pdb.name = 'in.pdb'
        t.get_current_artifact = lambda k: state
        seen = {}

        def model_one_tail(*a, **kw):
            seen['seed'] = kw['seed']
            raise self._Stop()
        with mock.patch('pestifer.psfutil.loop_ccd.loop_atoms_from_pdb', return_value=([], [], [])), \
             mock.patch('pestifer.psfutil.tail_model.model_one_tail', side_effect=model_one_tail):
            with self.assertRaises(self._Stop):
                t.model_terminal_tails()
        return seen['seed'], t.outcome

    def test_a_build_without_seed_keeps_the_old_seed(self):
        """Default builds must be unchanged: namd.seed at its default gives exactly 90125."""
        legacy = BaseTask._LEGACY_SEED
        self.assertEqual(self._task(legacy).declash_seed(), 90125)
        self.assertEqual(self._task(None).declash_seed(), 90125)

    def test_namd_seed_zero_does_not_go_negative(self):
        self.assertGreaterEqual(self._task(0).declash_seed(), 0)

    def test_the_tail_modeler_takes_the_build_seed(self):
        seed, outcome = self._tail_seed(27021973)
        self.assertEqual(seed, 90126, 'the first tail must model from the declash seed + 0')
        self.assertEqual(outcome.get('declash_seed'), 90126)
        self.assertEqual(self._tail_seed(27021972)[0], 90125)

    def test_the_pendant_declasher_takes_the_build_seed(self):
        expect = self._first_draw(np.random.default_rng(90126))
        self.assertEqual(self._pendant_rng(27021973), expect)
        self.assertNotEqual(self._pendant_rng(27021973), self._pendant_rng(27021974))

    def test_the_loop_declasher_takes_the_build_seed(self):
        expect = self._first_draw(np.random.default_rng(90126))
        self.assertEqual(self._loop_rng(27021973), expect)
        self.assertNotEqual(self._loop_rng(27021973), self._loop_rng(27021974))


if __name__ == '__main__':
    unittest.main()
