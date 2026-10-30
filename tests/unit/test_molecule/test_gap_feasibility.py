# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""An interior gap that cannot close against its deposited anchors is widened by releasing flanking
resolved residues into the loop (``Molecule.widen_unclosable_gaps``), and psfgen says so.

8FAD lists gp120 N141-N142-S143 (label B/D/E 109-111) as unobserved, but the deposited I140 and S144
point their peptide planes so that three residues cannot join them: built as deposited, Ser112 came
out cis and then D in all three chains of example 9.  Releasing 108 and 112 lets the five close trans
and L, and the rebuilt anchors sit within 0.9 A of their deposited positions.
"""
import inspect
import logging
import os
import re
import shutil
import tempfile
import unittest
from types import SimpleNamespace as NS

import numpy as np

from pestifer.molecule.molecule import Molecule
from pestifer.molecule.stateinterval import StateInterval
from pestifer.psfutil.loop_ccd import gap_closes, place_atom_nerf, _TRIAL
from pestifer.tasks.psfgen import PsfgenTask

HERE = os.path.dirname(os.path.abspath(__file__))


def ideal_chain(nres, phipsi=(-63.0, -43.0)):
    """Backbone dicts (N, CA, C, O) of ``nres`` ideal residues with one (phi, psi)."""
    t = _TRIAL
    N, CA = np.array([0., 0., 0.]), np.array([t['N_CA'], 0., 0.])
    C = CA + t['CA_C'] * np.array([-np.cos(np.radians(t['N_CA_C'])), np.sin(np.radians(t['N_CA_C'])), 0.])
    out = []
    for i in range(nres):
        nN = place_atom_nerf(N, CA, C, t['C_N'], t['CA_C_N'], phipsi[1])
        O = place_atom_nerf(nN, CA, C, 1.23, 120.5, 180.0)   # in the peptide plane, opposite N(i+1)
        out.append({'N': N, 'CA': CA, 'C': C, 'O': O})
        nCA = place_atom_nerf(CA, C, nN, t['N_CA'], t['C_N_CA'], 180.0)
        nC = place_atom_nerf(C, nN, nCA, t['CA_C'], t['N_CA_C'], phipsi[0])
        N, CA, C = nN, nCA, nC
    return out


class TestGapCloses(unittest.TestCase):

    def test_anchors_cut_from_a_real_chain_close(self):
        # POSITIVE CONTROL: a gap opened in an ideal helix can be bridged by the residues removed
        chain = ideal_chain(12)
        for n in (1, 3, 6):
            ok, rmsd = gap_closes(chain[2], chain[3 + n], n, np.random.default_rng(0))
            self.assertTrue(ok, f'{n} residues, best {rmsd:.2f} A')

    def test_anchors_out_of_reach_do_not(self):
        chain = ideal_chain(8)
        far = {k: v + np.array([40., 0., 0.]) for k, v in chain[6].items()}
        ok, rmsd = gap_closes(chain[2], far, 3, np.random.default_rng(0))
        self.assertFalse(ok)
        self.assertGreater(rmsd, 5.0)

    def test_the_n_anchor_carbonyl_is_not_free(self):
        """Turning the N-anchor's peptide plane away (its O, which fixes its psi) is what stops
        8FAD's gaps; flipping it here must cost closure that the untouched anchor allows."""
        chain = ideal_chain(8)
        A = dict(chain[2])
        ok, _ = gap_closes(A, chain[4], 1, np.random.default_rng(0))
        self.assertTrue(ok)
        axis = A['C'] - A['CA']
        axis /= np.linalg.norm(axis)
        v = A['O'] - A['C']
        A['O'] = A['C'] + (2 * np.dot(v, axis) * axis - v)   # rotate O 180 deg about CA-C
        self.assertFalse(gap_closes(A, chain[4], 1, np.random.default_rng(0))[0])


class TestBuildsAsLoop(unittest.TestCase):

    def test_length_or_force(self):
        short = StateInterval(state='MISSING', bounds=[3, 5])
        self.assertFalse(short.builds_as_loop(4))
        self.assertTrue(short.builds_as_loop(3))
        short.force_loop = True
        self.assertTrue(short.builds_as_loop(4), 'a widened gap is a loop however short')
        self.assertFalse(StateInterval(state='RESOLVED', bounds=[0, 9], force_loop=True).builds_as_loop(1))

    def test_no_build_path_compares_a_gap_length_directly(self):
        """Every place that decides loop-or-continuous must go through builds_as_loop, or a widened
        3-residue gap is built continuously at that site (the defect this exists to fix)."""
        root = os.path.dirname(os.path.dirname(inspect.getfile(Molecule)))
        offenders = []
        for d, _dirs, files in os.walk(root):
            for fn in files:
                if fn.endswith('.py') and fn != 'stateinterval.py':    # where builds_as_loop lives
                    p = os.path.join(d, fn)
                    for k, ln in enumerate(open(p), 1):
                        if re.search(r'num_items\(\)\s*>=\s*min_(loop_)?length', ln) and not ln.lstrip().startswith('#'):
                            offenders.append(f'{p}:{k}')
        self.assertEqual(offenders, [])


class TestRealStructures(unittest.TestCase):
    """8FAD (fetched) needs two residues released in each gp120 chain -- 107 and 108, the N side
    winning each step on B-factor (Asn107 268, Ile108 254, Ser112 253); 4ZMJ (a tracked
    fixture), whose gaps are ordinary loops, needs nothing.  Each test parses its own copy:
    widening mutates the molecule."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.addCleanup(os.chdir, os.getcwd())
        os.chdir(self.d)

    def _molecule(self, sid):
        from pestifer.molecule.chainidmanager import ChainIDManager
        if sid == '8fad' and not os.path.exists('8fad.cif'):
            from pestifer.core.pipeline import PipelineContext
            from pestifer.tasks.fetch import FetchTask
            task = FetchTask(specs={'source': 'rcsb', 'sourceID': '8fad', 'source_format': 'cif'})
            task.provision(packet=dict(pipeline=PipelineContext()))
            task.execute()
        if sid == '4zmj':
            shutil.copy(os.path.join(HERE, 'test_molecule', '4zmj.cif'), self.d)
        return Molecule(source={'source_id': sid, 'source_db': 'rcsb', 'file_format': 'mmCIF',
                                'sequence': {'fix_conflicts': True, 'fix_engineered_mutations': True,
                                             'include_terminal_loops': False}},
                        molid=0, chainIDmanager=ChainIDManager(format='mmCIF'))

    @staticmethod
    def _segment(m, chain):
        return next(S for S in m.asymmetric_unit.segments if S.segtype == 'protein' and S.chainID == chain)

    def test_4zmj_releases_nothing(self):
        m = self._molecule('4zmj')
        nint = sum(1 for S in m.asymmetric_unit.segments if S.segtype == 'protein'
                   for i, b in enumerate(S.subsegments)
                   if b.state == 'MISSING' and 0 < i < len(S.subsegments) - 1)
        self.assertGreaterEqual(nint, 3, 'POSITIVE CONTROL: there must be gaps to judge')
        self.assertEqual(m.widen_unclosable_gaps(), [])

    def test_8fad_widens_each_gp120_gap_by_two(self):
        m = self._molecule('8fad')
        before = {c: [b.bounds for b in self._segment(m, c).subsegments] for c in 'BDE'}
        found = m.widen_unclosable_gaps()
        self.assertEqual(sorted((f['segment'].chainID, [r.resid.resid for r in f['gap']],
                                 [r.resid.resid for r in f['released']], f['closes']) for f in found),
                         [(c, [109, 110, 111], [107, 108], True) for c in 'BDE'])
        for c in 'BDE':
            S = self._segment(m, c)
            self.assertNotEqual([b.bounds for b in S.subsegments], before[c])
            loop = [b for b in S.subsegments if b.state == 'MISSING'
                    and [r.resid.resid for r in S.residues[b.bounds[0]:b.bounds[1] + 1]] == [107, 108, 109, 110, 111]]
            self.assertEqual(len(loop), 1, 'one five-residue gap from 107 to 111')
            self.assertTrue(loop[0].force_loop)
            self.assertTrue(loop[0].builds_as_loop(4))
            released = [r for r in S.residues if r.resid.resid in (107, 108)]
            self.assertTrue(all(not r.resolved and len(r.atoms) == 0 for r in released))

    def test_report_only_changes_nothing(self):
        m = self._molecule('8fad')
        bounds = [b.bounds for b in self._segment(m, 'B').subsegments]
        found = m.widen_unclosable_gaps(release=False)
        self.assertEqual(len(found), 3)
        self.assertEqual([b.bounds for b in self._segment(m, 'B').subsegments], bounds)
        self.assertTrue(all(r.resolved for r in self._segment(m, 'B').residues if r.resid.resid in (108, 112)))

    def test_a_bonded_anchor_is_never_released(self):
        m = self._molecule('8fad')
        topol = m.asymmetric_unit.objmanager.setdefault('topol', {})
        fake = NS(chainID1='B', resid1=NS(resid=108), chainID2='B', resid2=NS(resid=112))
        topol['ssbonds'] = list(topol.get('ssbonds', [])) + [fake]
        found = {f['segment'].chainID: f for f in m.widen_unclosable_gaps()}
        self.assertFalse(found['B']['closes'], 'with both anchors bonded, nothing can be released')
        self.assertEqual(found['B']['released'], [])
        self.assertTrue(all(r.resolved for r in self._segment(m, 'B').residues if r.resid.resid in (108, 112)))
        self.assertTrue(found['D']['closes'], 'NEGATIVE CONTROL: the other chains are untouched by it')

    def test_the_cap_holds(self):
        m = self._molecule('8fad')
        found = m.widen_unclosable_gaps(max_per_side=0)
        self.assertEqual(len(found), 3)
        self.assertTrue(all(not f['closes'] and f['released'] == [] for f in found))


class TestFollowedBy(unittest.TestCase):

    def test_only_later_tasks_count(self):
        t = PsfgenTask.__new__(PsfgenTask)
        t.index = 1
        t.provisions = {'tasknames': ['ligate', 'psfgen', 'md', 'ligate']}
        self.assertTrue(t.followed_by('ligate'))
        t.provisions = {'tasknames': ['ligate', 'psfgen', 'md']}
        self.assertFalse(t.followed_by('ligate'), 'an earlier ligate does not close a later loop')


class TestTheWarning(unittest.TestCase):

    @staticmethod
    def _res(chain, resid, auth_chain=None, auth=None):
        return NS(chainID=chain, resid=NS(resid=resid), auth_asym_id=auth_chain, auth_seq_id=auth)

    def _run(self, found, widen=True, ligate=True):
        calls = []
        task = PsfgenTask.__new__(PsfgenTask)
        task.specs = {'source': {'sequence': {'loops': {'widen_unclosable': widen}}}}
        task.followed_by = lambda name: ligate and name == 'ligate'
        task.base_molecule = NS(WIDEN_MAX_PER_SIDE=3,
                                widen_unclosable_gaps=lambda release: calls.append(release) or found)
        logname = 'pestifer.tasks.psfgen'
        with self.assertLogs(logname, level='DEBUG') as cm:
            logging.getLogger(logname).debug('-')
            task._widen_unclosable_gaps()
        return calls, [r.getMessage() for r in cm.records if r.levelno == logging.WARNING]

    def _found(self, closes=True):
        gap = [self._res('B', i, 'C', i + 32) for i in (109, 110, 111)]
        rel = [self._res('B', 108, 'C', 140), self._res('B', 112, 'C', 144)] if closes else []
        return [dict(segment=None, gap=gap, released=rel, closes=closes)]

    def test_a_release_is_named_in_both_numberings(self):
        calls, warn = self._run(self._found())
        self.assertEqual(calls, [True])
        self.assertEqual(len(warn), 1)
        self.assertIn('B:109-111 (author C:141-143)', warn[0])
        self.assertIn('released B:108 (author C:140), B:112 (author C:144)', warn[0])

    def test_without_ligate_nothing_is_released(self):
        calls, warn = self._run(self._found(), ligate=False)
        self.assertEqual(calls, [False])
        self.assertIn('not released because no ligate task follows', warn[0])

    def test_switched_off_nothing_is_released(self):
        calls, warn = self._run(self._found(), widen=False)
        self.assertEqual(calls, [False])
        self.assertIn('not released because widen_unclosable is false', warn[0])

    def test_a_gap_that_cannot_close_says_so(self):
        _calls, warn = self._run(self._found(closes=False))
        self.assertIn('cannot close against its resolved anchors, even with up to 3', warn[0])

    def test_nothing_to_widen_says_nothing(self):
        _calls, warn = self._run([])
        self.assertEqual(warn, [])

    def test_the_build_calls_it_before_psfgen_runs(self):
        src = inspect.getsource(PsfgenTask.do)
        call, run = src.find('self._widen_unclosable_gaps()'), src.find('self.result = self.psfgen()')
        self.assertTrue(0 < call < run, 'the gaps must be widened before psfgen builds them')


if __name__ == '__main__':
    unittest.main()
