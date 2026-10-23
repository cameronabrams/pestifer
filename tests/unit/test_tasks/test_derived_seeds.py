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

import yaml

from pestifer.tasks.basetask import BaseTask
from pestifer.tasks.ligate import LigateTask
from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask


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


if __name__ == '__main__':
    unittest.main()
