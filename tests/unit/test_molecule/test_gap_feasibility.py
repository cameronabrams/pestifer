# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""An interior gap whose resolved anchors cannot hold the residues declared missing
(``Molecule.infeasible_protein_gaps``) and the psfgen warning that names the fix.

8FAD lists B/D/E 109-111 as unobserved, but Ile108 and Ser112 sit 3.89 A apart (CA-CA), the spacing
of a bonded pair.  psfgen squeezes three residues in anyway; Ser112 comes out cis and then D at the
first minimize.  Built with the three deleted, 108-112 closes trans and L.
"""
import inspect
import logging
import os
import shutil
import tempfile
import unittest
from types import SimpleNamespace as NS

import numpy as np

from pestifer.molecule.molecule import Molecule
from pestifer.tasks.psfgen import PsfgenTask

STEP = Molecule.CA_STEP


def residue(resid, ca=None):
    atoms = [] if ca is None else [NS(name='N', x=0., y=0., z=0.),
                                   NS(name='CA', x=float(ca[0]), y=float(ca[1]), z=float(ca[2]))]
    return NS(resid=NS(resid=resid), atoms=NS(data=atoms))


def chain(span, nmissing, chainID='A'):
    """Resolved 1-2, then ``nmissing`` missing residues, then resolved anchors whose CAs sit
    ``span`` A from residue 2's."""
    res = [residue(1, (-STEP, 0, 0)), residue(2, (0, 0, 0))]
    res += [residue(3 + k) for k in range(nmissing)]
    c = 3 + nmissing
    res += [residue(c, (span, 0, 0)), residue(c + 1, (span + STEP, 0, 0))]
    sub = lambda state, a, b: NS(state=state, bounds=(a, b), num_items=lambda: b - a + 1)
    subs = [sub('RESOLVED', 0, 1), sub('MISSING', 2, 1 + nmissing),
            sub('RESOLVED', 2 + nmissing, 3 + nmissing)]
    return NS(segtype='protein', chainID=chainID, segname=chainID, subsegments=subs, residues=res)


def molecule(*segments):
    m = Molecule.__new__(Molecule)
    m.asymmetric_unit = NS(segments=list(segments))
    return m


class TestTheRule(unittest.TestCase):

    def test_an_ordinary_gap_is_not_reported(self):
        # NEGATIVE CONTROL for both sides: 3 residues across 10 A is a normal loop
        self.assertEqual(molecule(chain(10.0, 3)).infeasible_protein_gaps(), [])

    def test_anchors_at_bonded_spacing_are_too_close(self):
        (g,) = molecule(chain(3.89, 3, 'B')).infeasible_protein_gaps()
        self.assertEqual((g['problem'], g['chainID'], g['first'], g['last'], g['nmissing']),
                         ('too_close', 'B', 3, 5, 3))
        self.assertEqual((g['n_anchor'], g['c_anchor']), (2, 6))
        self.assertAlmostEqual(g['ca_ca'], 3.89)

    def test_the_close_threshold_is_inclusive_and_sharp(self):
        t = Molecule.GAP_ANCHORS_BONDED
        self.assertEqual(len(molecule(chain(t, 1)).infeasible_protein_gaps()), 1)
        self.assertEqual(molecule(chain(t + 0.01, 1)).infeasible_protein_gaps(), [])

    def test_anchors_beyond_reach_are_too_far(self):
        n = 2
        self.assertEqual(molecule(chain(STEP * (n + 1), n)).infeasible_protein_gaps(), [],
                         'fully extended is the limit, not past it')
        (g,) = molecule(chain(STEP * (n + 1) + 0.1, n)).infeasible_protein_gaps()
        self.assertEqual(g['problem'], 'too_far')

    def test_a_short_gap_is_checked_whatever_min_loop_length_is(self):
        m = molecule(chain(3.89, 1))
        m.min_loop_length = 4      # psfgen builds a gap below this continuously; still a problem
        self.assertEqual(len(m.infeasible_protein_gaps()), 1)

    def test_only_interior_protein_gaps_count(self):
        s = chain(3.89, 3)
        # a terminal gap has no anchor on one side
        terminal = chain(3.89, 3)
        terminal.subsegments = terminal.subsegments[1:]
        glycan = chain(3.89, 3)
        glycan.segtype = 'glycan'
        self.assertEqual(len(molecule(s, terminal, glycan).infeasible_protein_gaps()), 1)

    def test_an_anchor_without_a_ca_is_skipped(self):
        s = chain(3.89, 3)
        s.residues[1].atoms.data = s.residues[1].atoms.data[:1]
        self.assertEqual(molecule(s).infeasible_protein_gaps(), [])


class TestTheWarning(unittest.TestCase):

    def _warn(self, *segments):
        task = PsfgenTask.__new__(PsfgenTask)
        task.base_molecule = molecule(*segments)
        logname = 'pestifer.tasks.psfgen'
        with self.assertLogs(logname, level='DEBUG') as cm:
            logging.getLogger(logname).debug('-')
            task._warn_infeasible_gaps()
        return [r.getMessage() for r in cm.records if r.levelno == logging.WARNING]

    def test_too_close_names_the_deletion_to_paste(self):
        warn = self._warn(chain(3.89, 3, 'B'), chain(3.89, 3, 'D'))
        self.assertEqual(len(warn), 3)
        self.assertIn('B:3-5: 3 residue(s) declared missing', warn[0])
        self.assertIn('deletions: [B:3-5, D:3-5]', warn[-1])

    def test_too_far_does_not_suggest_a_deletion(self):
        warn = self._warn(chain(20.0, 2, 'C'))
        self.assertEqual(len(warn), 1)
        self.assertIn('cannot span', warn[0])
        self.assertNotIn('deletions', warn[0])

    def test_a_feasible_build_says_nothing(self):
        self.assertEqual(self._warn(chain(10.0, 3)), [])

    def test_the_build_calls_it_before_psfgen_runs(self):
        src = inspect.getsource(PsfgenTask.do)
        call, run = src.find('self._warn_infeasible_gaps()'), src.find('self.result = self.psfgen()')
        self.assertTrue(0 < call < run, 'the check must run, and run before psfgen builds the loop')


class TestRealStructures(unittest.TestCase):
    """8FAD (fetched) carries exactly the three impossible gaps; 4ZMJ (a tracked fixture), whose
    three gaps per protomer are ordinary loops, carries none."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        self.cwd = os.getcwd()
        self.addCleanup(os.chdir, self.cwd)

    def _molecule(self, sid):
        from pestifer.molecule.chainidmanager import ChainIDManager
        os.chdir(self.d)
        return Molecule(source={'source_id': sid, 'source_db': 'rcsb', 'file_format': 'mmCIF',
                                'sequence': {'fix_conflicts': True, 'fix_engineered_mutations': True,
                                             'include_terminal_loops': False}},
                        molid=0, chainIDmanager=ChainIDManager(format='mmCIF'))

    def test_4zmj_has_no_infeasible_gap(self):
        shutil.copy(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'test_molecule', '4zmj.cif'),
                    self.d)
        m = self._molecule('4zmj')
        nint = sum(1 for S in m.asymmetric_unit.segments if S.segtype == 'protein'
                   for i, b in enumerate(S.subsegments)
                   if b.state == 'MISSING' and 0 < i < len(S.subsegments) - 1)
        self.assertGreaterEqual(nint, 3, 'POSITIVE CONTROL: there must be gaps to judge')
        self.assertEqual(m.infeasible_protein_gaps(), [])

    def test_8fad_reports_b_d_e_109_111(self):
        from pestifer.core.pipeline import PipelineContext
        from pestifer.tasks.fetch import FetchTask
        os.chdir(self.d)
        task = FetchTask(specs={'source': 'rcsb', 'sourceID': '8fad', 'source_format': 'cif'})
        task.provision(packet=dict(pipeline=PipelineContext()))
        task.execute()
        found = self._molecule('8fad').infeasible_protein_gaps()
        self.assertEqual(sorted((g['chainID'], g['first'], g['last'], g['problem']) for g in found),
                         [(c, 109, 111, 'too_close') for c in 'BDE'])
        for g in found:
            self.assertTrue(np.isclose(g['ca_ca'], 3.89, atol=0.01))


if __name__ == '__main__':
    unittest.main()
