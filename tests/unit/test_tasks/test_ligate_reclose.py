# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A loop that still pierces a ring after rotation is re-closed with the other loops present.

Example 14 replica 1 in the 3.29.8 sweep: the two symmetry copies of the 37-residue loop 657-693
were closed in parallel, each blind to the other, and closed into each other -- a bond of
A SER 678 through the ring of E TYR 675.  No rotation could clear it, and terminate refused the
structure 1.5 hours of MD later.

Pure logic (the closure and the ring engine are mocked), so deliberately not ``needs_tools``:
CI has to run it.
"""
import logging
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


class TestARotationThatOpensAJunctionIsReclosed(unittest.TestCase):
    """Example 11 (7TXD), after the psfgen-time closure landed: the piercing check cleared HSD D-53
    by rotating a backbone branch of loop D:32-52 at Ile48.  That swung the loop's free end 6.8 A
    off its anchor, connect bonded 52-53 as it lay, and the junction minimized to a cis bond
    (omega -29 deg) in the final structure -- filed by the tracker as inherited."""

    def _run(self, t, engine_results, opened):
        from pestifer.tasks.ligate import LigateTask
        with mock.patch('pestifer.psfutil.loop_ccd.loop_atoms_from_pdb', side_effect=_loop_atoms), \
             mock.patch('pestifer.psfutil.loop_ccd.close_one_loop', side_effect=_closed) as close, \
             mock.patch('pestifer.psfutil.loop_ccd.heavy_coords_of_serials',
                        side_effect=lambda _p, s: np.ones((len(s), 3))), \
             mock.patch('pestifer.util.coord.pdb_replace_coords'), \
             mock.patch('pestifer.psfutil.ring_resolve.resolve_model_built_piercings',
                        side_effect=engine_results), \
             mock.patch.object(LigateTask, 'open_junctions', side_effect=opened) as oj, \
             self.assertLogs('pestifer.tasks.loopclosure', level='INFO') as cm:
            t.resolve_loop_piercings()
            logging.getLogger('pestifer.tasks.loopclosure').info('end')   # a silent pass is a pass
        return close, oj, [r.getMessage() for r in cm.records if r.levelname == 'WARNING']

    def test_the_opened_loop_is_reclosed_and_checked_again(self):
        t = _task()
        close, oj, warns = self._run(t, [(None, []), (None, [])], [{1: (6.82, -28.9)}, {}])
        close.assert_called_once()
        self.assertEqual(close.call_args.kwargs['segname'], 'A')
        self.assertEqual(oj.call_count, 2, 'the re-closed structure must be checked again')
        self.assertEqual(warns, [])
        self.assertEqual(t.outcome['loop_junctions_open'], 0)

    def test_one_that_stays_open_is_named(self):
        from pestifer.tasks.ligate import _RECLOSE_ATTEMPTS
        t = _task()
        still = {1: (6.82, -28.9)}
        close, _oj, warns = self._run(t, [(None, [])] * (_RECLOSE_ATTEMPTS + 1),
                                      [still] * (_RECLOSE_ATTEMPTS + 1))
        self.assertEqual(close.call_count, _RECLOSE_ATTEMPTS)
        self.assertTrue(any('A:657-659 no longer meets its anchor A:660' in w for w in warns), warns)
        self.assertEqual(t.outcome['loop_junctions_open'], 1)

    def test_steer_is_not_held_to_a_ccd_junction(self):
        """NEGATIVE CONTROL: a steered loop's end never sat on an ideal junction."""
        t = _task(with_ccd_context=False)
        close, oj, warns = self._run(t, [(None, [])], [{1: (6.82, -28.9)}])
        oj.assert_not_called()
        close.assert_not_called()


class TestOpenJunctions(unittest.TestCase):
    """The detector, on a real structure: BPTI's 28-29 peptide bond as the junction of a loop
    24-28, intact, pulled apart, and turned cis."""

    def setUp(self):
        import os, tempfile, shutil
        from pathlib import Path
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.src = str(Path(__file__).parent / 'fixtures' / 'continuation_inputs' / 'my_6pti.pdb')
        self.pdb = os.path.join(self.d, 'x.pdb')
        self.gaps = [{'segname': 'A', 'loop_resids': [24, 25, 26, 27, 28], 'n_anchor_resid': 23,
                      'c_anchor_resid': 29}]
        t = _task()
        state = mock.Mock()
        state.pdb.name = self.pdb
        t.get_current_artifact = lambda k: state
        self.t = t

    def _write(self, move):
        rows = []
        for ln in open(self.src):
            if ln.startswith('ATOM') and ln[72:76].strip() == 'A':
                rows.append(ln)
        xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in rows])
        key = [(int(l[22:26]), l[12:16].strip()) for l in rows]
        xyz = move(xyz, key)
        with open(self.pdb, 'w') as f:
            for l, (x, y, z) in zip(rows, xyz):
                f.write(f'{l[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{l[54:]}')

    def test_an_intact_junction_is_closed(self):
        self._write(lambda xyz, key: xyz)
        self.assertEqual(self.t.open_junctions(self.gaps), {})

    def test_an_end_swung_off_its_anchor_is_open(self):
        def swing(xyz, key):
            loop = [i for i, (r, _n) in enumerate(key) if 24 <= r <= 28]
            xyz[loop] += [4.0, 3.0, 0.0]
            return xyz
        self._write(swing)
        d, _w = self.t.open_junctions(self.gaps)[0]
        self.assertGreater(d, 4.0)

    def test_a_cis_junction_is_open_though_bonded(self):
        from pestifer.util.coord import rotate_points_about_axis

        def turn(xyz, key):
            i = {k: j for j, k in enumerate(key)}
            c, n = xyz[i[(28, 'C')]], xyz[i[(29, 'N')]]
            down = [j for j, (r, nm) in enumerate(key) if r > 29 or (r == 29 and nm != 'N')]
            xyz[down] = rotate_points_about_axis(xyz[down], n, n - c, 170.0)
            return xyz
        self._write(turn)
        got = self.t.open_junctions(self.gaps)
        self.assertIn(0, got)
        self.assertLess(abs(got[0][1]), 30.0)


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
