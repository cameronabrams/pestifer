# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Backbone knot detection (``psfutil/knots.py``) and its place in the stereochemistry tracker.

The invariant is checked against knots whose type is known by construction (unknot, trefoil,
figure-eight, the (2,5) torus knot) and against a real CHARMM build (BPTI), which must read
unknotted.  A knotted BPTI is made by translating each residue rigidly so its alpha carbon lands
on a trefoil with straight tails: rigid moves keep every stereocenter, so the knot is the only
new thing in the frame.  Real knotted proteins -- YbeA (1NS5, trefoil, determinant 3) and UCH-L3
(1XD3, 5_2, determinant 7) -- were checked by hand when this was written; they are not fixtures.
"""
import logging
import os
import shutil
import tempfile
import unittest
from types import SimpleNamespace as NS

import numpy as np

from pestifer.psfutil.knots import close_to_infinity, kmt_reduce, knot_determinant
from pestifer.psfutil.stereo import StereoChecker
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


class FakeTask:
    """Just enough of a task for ``BaseTask.track_stereochemistry``."""
    track_stereochemistry = BaseTask.track_stereochemistry

    def __init__(self, store, pdb, taskname):
        self.store, self.taskname = store, taskname
        self.pipeline = store.setdefault('_pipeline', NS())
        self.state = NS(psf=NS(path=PSF), pdb=NS(path=pdb))

    def get_current_artifact_data(self, key):
        return self.store.get(key)

    def get_current_artifact(self, key):
        return self.state if key == 'state' else None

    def register(self, data, key):
        self.store[key] = data

T = np.linspace(0, 2 * np.pi, 240, endpoint=False)
UNKNOT = np.c_[np.cos(T), np.sin(T), 0 * T]
TREFOIL = np.c_[np.sin(T) + 2 * np.sin(2 * T), np.cos(T) - 2 * np.cos(2 * T), -np.sin(3 * T)]
FIGURE8 = np.c_[(2 + np.cos(2 * T)) * np.cos(3 * T), (2 + np.cos(2 * T)) * np.sin(3 * T), np.sin(4 * T)]
TORUS25 = np.c_[np.cos(2 * T) * (3 + np.cos(5 * T)), np.sin(2 * T) * (3 + np.cos(5 * T)), np.sin(5 * T)]


def opened(K, tail=10.0, n=8):
    """``K`` cut at its outermost point, with straight tails leaving radially."""
    k = int(np.argmax(np.linalg.norm(K[:, :2], axis=1)))
    P = np.roll(K, -k - 1, axis=0)[:-2]
    c = P.mean(axis=0)
    d0 = (P[0] - c) / np.linalg.norm(P[0] - c)
    d1 = (P[-1] - c) / np.linalg.norm(P[-1] - c)
    return np.vstack([P[0] + np.outer(np.linspace(tail, 0.3, n), d0), P,
                      P[-1] + np.outer(np.linspace(0.3, tail, n), d1)])


def open_det(P):
    return knot_determinant(kmt_reduce(close_to_infinity(P)))


class TestTheInvariant(unittest.TestCase):

    def test_closed_knots_of_known_type(self):
        for K, det in ((UNKNOT, 1), (TREFOIL, 3), (FIGURE8, 5), (TORUS25, 5)):
            for seed in (1, 2, 3):            # any generic projection gives the same answer
                self.assertEqual(knot_determinant(K, seed=seed), det)

    def test_reduction_keeps_the_knot_and_collapses_the_unknot(self):
        self.assertEqual(len(kmt_reduce(UNKNOT)), 3)
        for K, det in ((TREFOIL, 3), (FIGURE8, 5)):
            R = kmt_reduce(K)
            self.assertLess(len(R), 20)
            self.assertEqual(knot_determinant(R), det)

    def test_an_open_chain_is_closed_far_outside(self):
        self.assertEqual(open_det(opened(TREFOIL)), 3)
        self.assertEqual(open_det(opened(FIGURE8)), 5)
        self.assertEqual(open_det(opened(UNKNOT)), 1, 'a cut circle is not knotted')

    def test_a_mirror_image_is_the_same_knot(self):
        self.assertEqual(open_det(opened(TREFOIL * [1, 1, -1])), 3)


def knotted_bpti(checker, xyz):
    """``xyz`` with each residue translated so its CA sits on an opened trefoil."""
    atoms = checker.topol.atoms.data
    chain = checker._chains[0] if getattr(checker, '_chains', None) else None
    if chain is None:
        checker.knots(xyz)
        chain = checker._chains[0]
    curve = opened(TREFOIL, tail=6.0, n=6)
    idx = np.linspace(0, len(curve) - 1, len(chain)).round().astype(int)
    target = 12.0 * curve[idx]                  # ~4 A between consecutive alpha carbons
    out = xyz.copy()
    res_of = {}
    for i, a in enumerate(atoms):
        res_of.setdefault((a.segname, a.resid.resid), []).append(i)
    for k, ca in enumerate(chain):
        rows = res_of[(atoms[ca].segname, atoms[ca].resid.resid)]
        out[rows] += target[k] - xyz[ca]
    return out


class TestOnARealBuild(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.checker = StereoChecker(PSF)

    def setUp(self):
        self.xyz = self.checker.load_coords(PDB)

    def test_bpti_is_one_unknotted_chain(self):
        self.assertEqual(self.checker.knots(self.xyz), [])
        atoms = self.checker.topol.atoms.data
        every_ca = [i for i, x in enumerate(atoms) if x.atomname == 'CA' and x.segtype == 'protein']
        self.assertEqual(self.checker._chains, [every_ca], 'one chain, every alpha carbon, in order')

    def test_a_knotted_bpti_is_named(self):
        kn = self.checker.knots(knotted_bpti(self.checker, self.xyz))
        self.assertEqual(len(kn), 1)
        self.assertEqual((kn[0].kind, kn[0].determinant, kn[0].first, kn[0].last), ('knot', 3, '1', '57'))
        self.assertIn('tied in a trefoil', kn[0].describe())

    def test_tying_the_knot_makes_no_stereocenter_defect(self):
        """So the tracker tests below are about the knot and nothing else."""
        k = knotted_bpti(self.checker, self.xyz)
        self.assertEqual([d for d in self.checker.check_coords(k) if d.kind != 'cis'], [])


class TestTheTrackerReportsAKnot(unittest.TestCase):

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, ignore_errors=True)
        ch = StereoChecker(PSF)
        self.knotted = os.path.join(self.d, 'knotted.pdb')
        write_pdb(PDB, self.knotted, knotted_bpti(ch, ch.load_coords(PDB)))
        self.store = {}

    def _run(self, pdb, name, **kw):
        with self.assertLogs('pestifer.tasks.basetask', level='DEBUG') as cm:
            logging.getLogger('pestifer.tasks.basetask').debug('-')
            FakeTask(self.store, pdb, name).track_stereochemistry(**kw)
        return [r for r in cm.records if r.getMessage() != '-']

    @staticmethod
    def _knot_warnings(recs):
        return [r.getMessage() for r in recs if r.levelno == logging.WARNING and 'trefoil' in r.getMessage()]

    def test_a_knot_is_blamed_on_the_stage_that_tied_it(self):
        self.assertEqual(self._knot_warnings(self._run(PDB, 'psfgen', baseline=True)), [])
        recs = self._run(self.knotted, 'md')
        self.assertTrue(self._knot_warnings(recs), [r.getMessage() for r in recs])
        self.assertTrue(any('md: introduced' in r.getMessage() for r in recs))
        recs = self._run(self.knotted, 'terminate', final=True)
        self.assertTrue(any('first seen at md' in m for m in self._knot_warnings(recs)))

    def test_a_knot_present_at_psfgen_still_warns(self):
        """Unlike a deposited cis bond, a knot as psfgen built it is not excused: the baseline
        already holds the loops psfgen placed, and proteins are almost never knotted."""
        self.assertTrue(self._knot_warnings(self._run(self.knotted, 'psfgen', baseline=True)))
        recs = self._run(self.knotted, 'terminate', final=True)
        self.assertTrue(any('still knotted as psfgen built it' in m for m in self._knot_warnings(recs)))

    def test_an_unknotted_build_says_nothing_about_knots(self):
        """NEGATIVE CONTROL for the two above: the same path on the real frame is quiet."""
        self._run(PDB, 'psfgen', baseline=True)
        recs = self._run(PDB, 'md') + self._run(PDB, 'terminate', final=True)
        self.assertEqual([r for r in recs if r.levelno >= logging.WARNING], [])


if __name__ == '__main__':
    unittest.main()
