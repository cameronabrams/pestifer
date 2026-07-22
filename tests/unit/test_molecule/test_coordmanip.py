import unittest
from pathlib import Path

import numpy as np

from pestifer.molecule.coordmanip import CoordManipulator, CoordManipulateError
from pestifer.objs.rottrans import RotTrans
from pestifer.objs.align import Align
from pestifer.objs.transfer_coords import TransferCoords

FIX = Path(__file__).parents[1] / 'test_tasks' / 'fixtures' / 'continuation_inputs'
PSF = str(FIX / 'my_6pti.psf')
PDB = str(FIX / 'my_6pti.pdb')

# atom counts confirmed against VMD atomselect on this exact fixture
_VMD_COUNTS = {
    'all': 14127, 'protein': 882, 'backbone': 229, 'name CA': 57,
    'chain A': 882, 'resid 10 to 20': 193, 'resid 42': 24,
    'protein and resid 5': 10, 'fragment 0': 882,
    'chain A and resid 1 to 5 and name CA': 5,
}


class TestSelection(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cm = CoordManipulator(PSF, PDB)

    def test_selection_counts_match_vmd(self):
        for sel, n in _VMD_COUNTS.items():
            self.assertEqual(int(self.cm.select(sel).sum()), n, f'selection {sel!r}')

    def test_and_is_intersection(self):
        a = self.cm.select('chain A')
        b = self.cm.select('name CA')
        self.assertTrue(np.array_equal(self.cm.select('chain A and name CA'), a & b))

    def test_resid_range_inclusive(self):
        self.assertEqual(int(self.cm.select('resid 5 to 5').sum()),
                         int(self.cm.select('resid 5').sum()))

    def test_unsupported_term_raises(self):
        for bad in ('resname ALA', 'not protein', 'chain A or chain B', 'within 5 of protein'):
            with self.assertRaises(CoordManipulateError):
                self.cm.select(bad)

    def test_fragment_and_disconnection_need_psf(self):
        ref = CoordManipulator(None, PDB)          # pdb-only: no bonds
        with self.assertRaises(CoordManipulateError):
            ref.select('fragment 0')


class TestRotTrans(unittest.TestCase):
    def setUp(self):
        self.cm = CoordManipulator(PSF, PDB)
        self.orig = self.cm.coords.copy()

    def test_trans_shifts_selection_only(self):
        # 'protein' is a fragment fully disconnected from the surrounding waters
        mask = self.cm.select('protein')
        self.cm.apply_rottrans(RotTrans(movetype='TRANS', x=10.0, y=-3.0, z=2.5, sel='protein'))
        d = self.cm.coords - self.orig
        self.assertTrue(np.allclose(d[mask], [10.0, -3.0, 2.5]))
        self.assertTrue(np.allclose(d[~mask], 0.0))

    def test_rot_is_rigid(self):
        # a rotation preserves all pairwise distances and the radius of gyration
        before = self.orig
        self.cm.apply_rottrans(RotTrans(movetype='ROT', axis='z', angle=57.0))
        after = self.cm.coords
        di = np.linalg.norm(before[:100] - before[0], axis=1)
        df = np.linalg.norm(after[:100] - after[0], axis=1)
        self.assertTrue(np.allclose(di, df, atol=1e-6))

    def test_rot_360_is_identity(self):
        self.cm.apply_rottrans(RotTrans(movetype='ROT', axis='y', angle=360.0))
        self.assertTrue(np.allclose(self.cm.coords, self.orig, atol=1e-6))

    def test_align_vector_onto_axis(self):
        # after ALIGN, the source direction is parallel to the target axis
        self.cm.apply_rottrans(RotTrans(movetype='ALIGN',
                                        source=['resid 5 and name CA', 'resid 50 and name CA'],
                                        target=[0.0, 0.0, 1.0]))
        a = self.cm.coords[self.cm.select('resid 5 and name CA')][0]
        b = self.cm.coords[self.cm.select('resid 50 and name CA')][0]
        v = (b - a) / np.linalg.norm(b - a)
        self.assertTrue(np.allclose(v, [0.0, 0.0, 1.0], atol=1e-5))

    def test_disconnection_guard(self):
        # a mid-chain residue range is bonded to its neighbors -> not a rigid fragment
        with self.assertRaises(CoordManipulateError):
            self.cm.apply_rottrans(RotTrans(movetype='TRANS', x=1.0, y=0.0, z=0.0,
                                            sel='resid 10 to 20'))


class TestAlignAndTransfer(unittest.TestCase):
    def test_align_onto_self_is_noop(self):
        cm = CoordManipulator(PSF, PDB)
        orig = cm.coords.copy()
        ref = CoordManipulator(PSF, PDB)
        cm.apply_align(Align(ref_pdb=PDB, ref_psf=PSF, mobile_sel='name CA',
                             ref_sel='name CA', apply_to='all'), ref)
        self.assertTrue(np.allclose(cm.coords, orig, atol=1e-4))

    def test_align_count_mismatch_raises(self):
        cm = CoordManipulator(PSF, PDB)
        ref = CoordManipulator(PSF, PDB)
        with self.assertRaises(CoordManipulateError):
            cm.apply_align(Align(ref_pdb=PDB, ref_psf=PSF, mobile_sel='name CA',
                                 ref_sel='resid 1 to 10 and name CA', apply_to='all'), ref)

    def test_transfer_copies_coords(self):
        cm = CoordManipulator(PSF, PDB)
        donor = CoordManipulator(PSF, PDB)
        donor.apply_rottrans(RotTrans(movetype='TRANS', x=5.0, y=5.0, z=5.0))  # move donor
        cm.apply_transfer_coords(TransferCoords(donor_pdb=PDB, donor_psf=PSF,
                                                donor_sel='resid 10 to 20',
                                                mobile_sel='resid 10 to 20'), donor)
        m = cm.select('resid 10 to 20')
        self.assertTrue(np.allclose(cm.coords[m], donor.coords[m]))

    def test_transfer_count_mismatch_raises(self):
        cm = CoordManipulator(PSF, PDB)
        donor = CoordManipulator(PSF, PDB)
        with self.assertRaises(CoordManipulateError):
            cm.apply_transfer_coords(TransferCoords(donor_pdb=PDB, donor_psf=PSF,
                                                    donor_sel='resid 10 to 20',
                                                    mobile_sel='resid 10 to 21'), donor)


if __name__ == '__main__':
    unittest.main()
