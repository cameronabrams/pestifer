# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Protein stereochemistry check (``psfutil/stereo.py``) and the build-stage tracker that uses it.

The oracle is a real CHARMM build (BPTI, ``my_6pti``), not a hand-made peptide: it must read
clean, which pins the sign conventions, and each defect is then made by moving its coordinates
the way a defect actually looks -- a mirror image for chirality, a rotation about one C-N bond for
a cis peptide.  Every defect test also asserts what must NOT be flagged, so a checker that flags
everything cannot pass.
"""
import logging
import os
import shutil
import tempfile
import unittest
from types import SimpleNamespace as NS

import numpy as np

from pestifer.psfutil.stereo import StereoChecker, CIS_CUTOFF
from pestifer.tasks.basetask import BaseTask

FIX = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   'test_tasks', 'fixtures', 'continuation_inputs')
PSF = os.path.join(FIX, 'my_6pti.psf')
PDB = os.path.join(FIX, 'my_6pti.pdb')


def write_pdb(src, dst, xyz):
    """Copy ``src`` with its ATOM/HETATM coordinates replaced, in file order."""
    i = 0
    with open(src) as f, open(dst, 'w') as g:
        for ln in f:
            if ln.startswith(('ATOM', 'HETATM')):
                x, y, z = xyz[i]
                ln = f'{ln[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{ln[54:]}'
                i += 1
            g.write(ln)


def rotate(xyz, rows, origin, axis, deg):
    axis = axis / np.linalg.norm(axis)
    t = np.radians(deg)
    p = xyz[rows] - origin
    xyz[rows] = (origin + p * np.cos(t) + np.cross(axis, p) * np.sin(t)
                 + np.outer(p @ axis, axis) * (1 - np.cos(t)))


class TestStereoChecker(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.checker = StereoChecker(PSF)
        cls.xyz = cls.checker.load_coords(PDB)

    def setUp(self):
        self.xyz = self.xyz.copy()          # each test perturbs its own copy

    def test_a_real_build_reads_clean(self):
        c = self.checker
        # POSITIVE CONTROL: BPTI has 58 residues, 6 Gly, 2 Thr... a checker that found no centers
        # or no peptide bonds would read clean for the wrong reason
        self.assertGreater(sum(1 for t in c.centers if t[0] == 'D-CA'), 45)
        self.assertGreater(sum(1 for t in c.centers if t[0] in ('THR-CB', 'ILE-CB')), 3)
        self.assertGreater(len(c.peptides), 50)
        self.assertEqual(c.check_coords(self.xyz), [])

    def test_a_mirror_image_is_all_D(self):
        self.xyz[:, 0] *= -1
        kinds = [d.kind for d in self.checker.check_coords(self.xyz)]
        ncenters = len(self.checker.centers)
        self.assertEqual(len([k for k in kinds if k != 'cis']), ncenters,
                         'every stereocenter inverts in a mirror image')
        self.assertNotIn('cis', kinds, 'a mirror image keeps every peptide trans')

    def _peptide(self, pro):
        """A peptide bond whose second residue is (or is not) proline."""
        for t in self.checker.peptides:
            if (t[2] == 'PRO') == pro:
                return t
        self.skipTest('no such peptide bond in the fixture')

    def _make_cis(self, t):
        """Rotate everything after C(i)-N(j), within the segment, so omega becomes ~0."""
        seg, resid, resname, ca_i, c_i, n_j, ca_j = t
        rows = [r for r in range(n_j + 1, self.checker.natoms)
                if self.checker.topol.atoms.data[r].segname == seg]
        from pestifer.psfutil.stereo import _dihedral
        w = float(_dihedral(*(self.xyz[[ca_i]], self.xyz[[c_i]], self.xyz[[n_j]],
                              self.xyz[[ca_j]]))[0])
        rotate(self.xyz, rows, self.xyz[n_j], self.xyz[n_j] - self.xyz[c_i], -w)
        return resid

    def test_a_cis_peptide_is_reported_at_its_residue(self):
        t = self._peptide(pro=False)
        resid = self._make_cis(t)
        found = self.checker.check_coords(self.xyz)
        self.assertEqual([(d.kind, d.resid) for d in found], [('cis', resid)],
                         'exactly the rotated bond, and no chirality change from a rigid rotation')
        self.assertLess(abs(found[0].value), CIS_CUTOFF)

    def test_a_cis_bond_before_proline_is_not_a_defect(self):
        self._make_cis(self._peptide(pro=True))
        self.assertEqual(self.checker.check_coords(self.xyz), [])

    def test_one_inverted_alpha_carbon(self):
        """Reflect one residue's CB through its N-CA-C plane: one D, nothing else."""
        kind, seg, resid, resname, ca, (n, c, cb) = next(
            t for t in self.checker.centers if t[0] == 'D-CA' and t[3] not in ('THR', 'ILE'))
        nrm = np.cross(self.xyz[n] - self.xyz[ca], self.xyz[c] - self.xyz[ca])
        nrm /= np.linalg.norm(nrm)
        self.xyz[cb] -= 2 * np.dot(self.xyz[cb] - self.xyz[ca], nrm) * nrm
        self.assertEqual([(d.kind, d.resid) for d in self.checker.check_coords(self.xyz)],
                         [('D-CA', resid)])

    def test_an_alpha_hydrogen_on_the_beta_side(self):
        """Reflect one HA through its N-CA-C plane: CB is untouched, so the residue is still L, and
        the only report is the hydrogen -- the state a coordinate transfer left 25-30 residues in."""
        seg, resid, resname, ca, n, c, cb, ha = self.checker.alpha_h[7]
        nrm = np.cross(self.xyz[n] - self.xyz[ca], self.xyz[c] - self.xyz[ca])
        nrm /= np.linalg.norm(nrm)
        self.xyz[ha] -= 2 * np.dot(self.xyz[ha] - self.xyz[ca], nrm) * nrm
        self.assertEqual([(d.kind, d.resid) for d in self.checker.check_coords(self.xyz)],
                         [('HA-CB', resid)])

    def test_a_real_d_residue_is_not_also_an_ha_defect(self):
        """NEGATIVE CONTROL for the above: a whole residue mirrored keeps HA opposite CB."""
        seg, resid, resname, ca, n, c, cb, ha = self.checker.alpha_h[7]
        nrm = np.cross(self.xyz[n] - self.xyz[ca], self.xyz[c] - self.xyz[ca])
        nrm /= np.linalg.norm(nrm)
        for i in (cb, ha):
            self.xyz[i] -= 2 * np.dot(self.xyz[i] - self.xyz[ca], nrm) * nrm
        kinds = {d.kind for d in self.checker.check_coords(self.xyz) if d.resid == resid}
        self.assertNotIn('HA-CB', kinds)

    def test_the_file_reader_and_the_array_agree(self):
        with tempfile.TemporaryDirectory() as d:
            self.xyz[:, 0] *= -1
            p = os.path.join(d, 'mirror.pdb')
            write_pdb(PDB, p, self.xyz)
            self.assertEqual(len(self.checker.check(p)), len(self.checker.centers))


class FakeTask:
    """Just enough of a task for ``BaseTask.track_stereochemistry``."""
    track_stereochemistry = BaseTask.track_stereochemistry

    def __init__(self, store, pdb, taskname, basename=None):
        self.store, self.taskname = store, taskname
        if basename:
            self.basename = basename
        self.pipeline = store.setdefault('_pipeline', NS())
        self.state = NS(psf=NS(path=PSF), pdb=NS(path=pdb))

    def get_current_artifact_data(self, key):
        return self.store.get(key)

    def get_current_artifact(self, key):
        return self.state if key == 'state' else None

    def register(self, data, key):
        self.store[key] = data


class TestTrackingAcrossStages(unittest.TestCase):
    """The tracker names the stage that introduced a defect, once, and only then."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        xyz = StereoChecker(PSF).load_coords(PDB)
        xyz[:, 0] *= -1
        self.mirror = os.path.join(self.d, 'mirror.pdb')
        write_pdb(PDB, self.mirror, xyz)
        self.store = {}

    def _run(self, pdb, name, **kw):
        with self.assertLogs('pestifer.tasks.basetask', level='DEBUG') as cm:
            logging.getLogger('pestifer.tasks.basetask').debug('-')
            FakeTask(self.store, pdb, name).track_stereochemistry(**kw)
        return [r for r in cm.records if r.getMessage() != '-']

    def test_a_defect_is_blamed_on_the_stage_that_made_it(self):
        recs = self._run(PDB, 'psfgen', baseline=True)
        self.assertFalse([r for r in recs if r.levelno >= logging.WARNING])
        self.assertEqual(self.store['stereo_record'], {})
        recs = self._run(self.mirror, 'md')
        warn = [r.getMessage() for r in recs if r.levelno == logging.WARNING]
        self.assertTrue(warn and 'md: introduced' in warn[0], warn)
        # reported once: the next stage carries it but did not make it
        recs = self._run(self.mirror, 'md2')
        self.assertFalse([r for r in recs if r.levelno >= logging.WARNING])
        recs = self._run(self.mirror, 'terminate', final=True)
        warn = [r.getMessage() for r in recs if r.levelno == logging.WARNING]
        self.assertTrue(any('first seen at md' in w for w in warn), warn)

    def test_a_defect_present_at_psfgen_is_not_blamed_on_anything(self):
        recs = self._run(self.mirror, 'psfgen', baseline=True)
        self.assertFalse([r for r in recs if r.levelno >= logging.WARNING],
                         'what psfgen hands on is inherited or placed, not introduced')
        self.assertTrue(any('present as psfgen wrote' in r.getMessage() for r in recs))
        recs = self._run(self.mirror, 'terminate', final=True)
        self.assertFalse([r for r in recs if r.levelno >= logging.WARNING])

    def test_the_stage_is_named_not_just_its_task(self):
        """A build runs several md tasks; "md: introduced" does not say which."""
        self._run(PDB, 'psfgen', baseline=True)
        with self.assertLogs('pestifer.tasks.basetask', level='WARNING') as cm:
            FakeTask(self.store, self.mirror, 'md', basename='00-04-000_md-minimize').track_stereochemistry()
        self.assertTrue(cm.output[0].endswith('md (00-04-000_md-minimize): introduced '
                                              f'{len(cm.output) - 1} protein stereochemistry defect(s) '
                                              'that no build step asked for:'), cm.output[0])
        recs = self._run(self.mirror, 'terminate', final=True)
        self.assertTrue(any('first seen at 00-04-000_md-minimize' in r.getMessage() for r in recs))

    def test_an_inherited_defect_that_comes_back_is_not_blamed(self):
        """A deposited omega near the cutoff wobbles across it in MD (4ZMJ G Thr139 at -24.5 deg):
        gone at one stage, back at the next.  Its return is not that stage's doing."""
        self._run(self.mirror, 'psfgen', baseline=True)
        self._run(PDB, 'md')                                  # everything inherited goes away
        recs = self._run(self.mirror, 'md2')                  # ... and comes back
        self.assertFalse([r for r in recs if r.levelno >= logging.WARNING],
                         [r.getMessage() for r in recs if r.levelno >= logging.WARNING][:2])
        self.assertTrue(any('again carries' in r.getMessage() for r in recs))
        recs = self._run(self.mirror, 'terminate', final=True)
        self.assertFalse([r for r in recs if r.levelno >= logging.WARNING])

    def test_an_inherited_defect_gone_at_the_end_is_not_counted(self):
        self._run(self.mirror, 'psfgen', baseline=True)
        self._run(PDB, 'md')
        recs = self._run(PDB, 'terminate', final=True)
        self.assertFalse([r for r in recs if r.levelno >= logging.WARNING])
        self.assertTrue(any('no protein stereochemistry defect was introduced' in r.getMessage()
                            for r in recs))

    def test_no_baseline_means_no_mid_build_check(self):
        """NEGATIVE CONTROL: without psfgen in this controller (a membrane patch's md stages),
        an md stage must not even read the structure -- the mirror would otherwise warn."""
        recs = self._run(self.mirror, 'md')
        self.assertEqual(recs, [])
        self.assertNotIn('stereo_record', self.store)

    def test_a_check_that_cannot_read_never_raises(self):
        self._run(PDB, 'psfgen', baseline=True)
        recs = self._run(os.path.join(self.d, 'missing.pdb'), 'md')
        self.assertTrue(any('skipped' in r.getMessage() for r in recs))


if __name__ == '__main__':
    unittest.main()


class TestTheBuildCallsIt(unittest.TestCase):
    """Aimed at the call sites: a tracker nothing calls is the original bug."""

    def test_an_md_task_checks_after_a_successful_run(self):
        from pestifer.tasks.mdtask import MDTask
        calls = []
        t = MDTask.__new__(MDTask)
        t.namdrun = lambda: 0
        t.track_stereochemistry = lambda **kw: calls.append(kw)
        self.assertEqual(t.do(), 0)
        self.assertEqual(calls, [{}])

    def test_a_failed_md_run_is_not_checked(self):
        from pestifer.tasks.mdtask import MDTask
        calls = []
        t = MDTask.__new__(MDTask)
        t.namdrun = lambda: 1
        t.track_stereochemistry = lambda **kw: calls.append(kw)
        t.do()
        self.assertEqual(calls, [])

    def test_a_loop_closure_checks_the_junctions_it_just_made(self):
        """Whichever task closes the loops (psfgen or ligate), the check comes after connect,
        when the junctions first exist -- and the molecule is marked closed only on success."""
        from pestifer.tasks.ligate import LigateTask
        calls = []
        t = LigateTask.__new__(LigateTask)
        t.basename = '00-01-004_psfgen-heal'
        t.base_molecule = NS(loops_closed_by=None)
        t.resolve_loop_piercings = lambda: calls.append('resolve')
        t.connect = lambda specs: calls.append('connect') or 0
        t.track_stereochemistry = lambda **kw: calls.append(('track', t.base_molecule.loops_closed_by))
        self.assertEqual(t.finish_loop_closure({}), 0)
        self.assertEqual(calls, ['resolve', 'connect', ('track', '00-01-004_psfgen-heal')])
        calls.clear()
        t.base_molecule.loops_closed_by = None
        t.connect = lambda specs: calls.append('connect') or 1
        self.assertEqual(t.finish_loop_closure({}), 1)
        self.assertIsNone(t.base_molecule.loops_closed_by, 'a failed connect closes nothing')

    def test_psfgen_sets_the_baseline_and_terminate_summarizes(self):
        import inspect
        from pestifer.tasks.psfgen import PsfgenTask
        from pestifer.tasks.terminate import TerminateTask
        self.assertEqual(inspect.getsource(PsfgenTask.do).count(
            'self.track_stereochemistry(baseline=True)'), 2,
            'both psfgen routes (build and readpsf-preserve) must record the baseline')
        self.assertIn('self.track_stereochemistry(final=True)', inspect.getsource(TerminateTask.do))
