# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A loop that still pierces a ring after rotation is re-closed with the other loops present.

Example 14 replica 1 in the 3.29.8 sweep: the two symmetry copies of the 37-residue loop 657-693
were closed in parallel, each blind to the other, and closed into each other -- a bond of
A SER 678 through the ring of E TYR 675.  No rotation could clear it, and terminate refused the
structure 1.5 hours of MD later.

Pure logic (the closure and the ring engine are mocked), so deliberately not ``needs_tools``:
CI has to run it.
"""
import unittest
from unittest import mock

import numpy as np

from pestifer.psfutil.loop_ccd import heavy_coords_of_serials

_GAPS = [
    {'segname': 'E', 'loop_resids': [657, 658, 659], 'n_anchor_resid': 656, 'c_anchor_resid': 660},
    {'segname': 'A', 'loop_resids': [657, 658, 659], 'n_anchor_resid': 656, 'c_anchor_resid': 660},
]
_SERIALS = {'E': [1, 2, 3], 'A': [4, 5, 6]}
_PIERCING = {'piercee': {'segname': 'E', 'resid': 658, 'resname': 'TYR'},
             'piercer': {'segname': 'A', 'resid': 659}}


def _task(with_ccd_context=True):
    from pestifer.tasks.ligate import LigateTask
    t = LigateTask.__new__(LigateTask)
    t.taskname, t.basename, t.specs, t.outcome = 'ligate', 'test', {}, {}
    t.next_basename = mock.Mock()
    t.register = mock.Mock()
    t.pdb_to_coor = mock.Mock()
    state = mock.Mock()
    state.psf.name, state.pdb.name = 'in.psf', 'in.pdb'
    state.xsc = None
    t.get_current_artifact = lambda k: state
    t.base_molecule = mock.Mock()
    t.base_molecule.protein_loop_gaps.return_value = _GAPS
    t.get_current_artifact_data = lambda k: t.base_molecule
    if with_ccd_context:
        t._ccd_context = dict(gaps=_GAPS, all_loop_serials=[1, 2, 3, 4, 5, 6], seed=100,
                              kwargs=dict(ensemble=2, refine=0, guard=False, max_iters=10, tol=0.1))
    return t


def _loop_atoms(_pdb, resids, segname=None):
    return None, None, _SERIALS[segname]


def _closed(**kw):
    ser = _SERIALS[kw['segname']]
    return {'serials': ser, 'closed': np.zeros((len(ser), 3))}


class TestReclose(unittest.TestCase):

    def _run(self, t, engine_results):
        env_calls = []

        def env(_pdb, serials):
            env_calls.append(sorted(serials))
            return np.ones((len(serials), 3))

        with mock.patch('pestifer.psfutil.loop_ccd.loop_atoms_from_pdb', side_effect=_loop_atoms), \
             mock.patch('pestifer.psfutil.loop_ccd.close_one_loop', side_effect=_closed) as close, \
             mock.patch('pestifer.psfutil.loop_ccd.heavy_coords_of_serials', side_effect=env), \
             mock.patch('pestifer.util.coord.pdb_replace_coords'), \
             mock.patch('pestifer.psfutil.ring_resolve.resolve_model_built_piercings',
                        side_effect=engine_results) as engine:
            t.resolve_loop_piercings()
        return close, engine, env_calls

    def test_the_piercing_loop_is_reclosed_against_the_other_loop(self):
        t = _task()
        # rotation fails, then the re-closed structure is clean
        close, engine, env_calls = self._run(t, [(None, [_PIERCING]), (None, [])])
        close.assert_called_once()
        kw = close.call_args.kwargs
        self.assertEqual(kw['segname'], 'A', 'the loop carrying the piercing bond is re-closed')
        self.assertEqual(env_calls, [_SERIALS['E']],
                         'the other loop must be in the environment: that is the whole fix')
        self.assertEqual(len(kw['extra_env']), 3)
        self.assertNotIn(kw['seed'], (100, 101), 'a re-closure must not redraw the first pose')
        self.assertEqual(engine.call_count, 2, 'the re-closed structure must be checked again')
        self.assertEqual(t.outcome['loop_piercings_unresolved'], 0)
        self.assertEqual(t.outcome['loop_reclosures'][0]['loop'], 'A:657-659')

    def test_the_reclosed_structure_becomes_the_state(self):
        t = _task()
        self._run(t, [(None, [_PIERCING]), (None, [])])
        self.assertTrue(any(c.kwargs.get('key') == 'state' for c in t.register.call_args_list))

    def test_retries_from_fresh_seeds_then_gives_up(self):
        from pestifer.tasks.ligate import _RECLOSE_ATTEMPTS
        t = _task()
        results = [(None, [_PIERCING])] * (_RECLOSE_ATTEMPTS + 1)
        close, engine, _env = self._run(t, results)
        seeds = [c.kwargs['seed'] for c in close.call_args_list]
        self.assertEqual(len(seeds), _RECLOSE_ATTEMPTS)
        self.assertEqual(len(set(seeds)), _RECLOSE_ATTEMPTS, 'each attempt needs its own seed')
        self.assertEqual(t.outcome['loop_piercings_unresolved'], 1)

    def test_a_ring_in_a_loop_pierced_from_outside_recloses_the_ring_loop(self):
        t = _task()
        outside = {'piercee': {'segname': 'E', 'resid': 658}, 'piercer': {'segname': 'B', 'resid': 9}}
        close, _e, _env = self._run(t, [(None, [outside]), (None, [])])
        self.assertEqual(close.call_args.kwargs['segname'], 'E')

    def test_steer_has_nothing_to_reclose(self):
        """Negative control: only a ccd closure leaves the context a re-closure needs."""
        t = _task(with_ccd_context=False)
        close, engine, _env = self._run(t, [(None, [_PIERCING])])
        close.assert_not_called()
        self.assertEqual(engine.call_count, 1)

    def test_a_clean_closure_recloses_nothing(self):
        t = _task()
        close, engine, _env = self._run(t, [(None, [])])
        close.assert_not_called()


class TestHeavyCoordsOfSerials(unittest.TestCase):

    def test_selects_exactly_the_named_heavy_atoms(self):
        import tempfile, os
        lines = []
        for serial, name, x in [(1, 'N', 1.0), (2, 'HN', 2.0), (3, 'CA', 3.0), (4, 'C', 4.0)]:
            lines.append(f'ATOM  {serial:>5d} {name:<4s} ALA A   1    {x:8.3f}{0.0:8.3f}{0.0:8.3f}'
                         f'  1.00  0.00      A\n')
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, 'x.pdb')
            with open(p, 'w') as f:
                f.writelines(lines)
            got = heavy_coords_of_serials(p, [2, 3, 4])
        self.assertEqual(got[:, 0].tolist(), [3.0, 4.0])


if __name__ == '__main__':
    unittest.main()
