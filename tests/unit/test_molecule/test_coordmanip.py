import unittest
from pathlib import Path

import numpy as np

from pestifer.molecule.coordmanip import CoordManipulator, CoordManipulateError
from pestifer.objs.rottrans import RotTrans
from pestifer.objs.align import Align
from pestifer.objs.transfer_coords import TransferCoords
from pestifer.objs.crot import Crot
from pestifer.objs.resid import ResID
from pestifer.psfutil.loop_ccd import dihedral_deg

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


def _phi(cm, chain, resid):
    """Measure the phi dihedral (C(-1)-N-CA-C) of chain/resid from a CoordManipulator."""
    cm._residue_index()
    r = cm._residue_of(chain, resid)
    coords, ridx, names = cm.coords, cm._ridx, cm._names
    p = [cm._res_xyz(coords, ridx, names, rr, nm)
         for rr, nm in ((r - 1, 'C'), (r, 'N'), (r, 'CA'), (r, 'C'))]
    return dihedral_deg(*p)


class TestCrot(unittest.TestCase):
    def setUp(self):
        self.cm = CoordManipulator(PSF, PDB)
        self.orig = self.cm.coords.copy()

    def test_residue_index_matches_resid_offset(self):
        # 6PTI is one chain numbered from 1, so VMD residue index == resid - 1
        self.assertEqual(self.cm._residue_of('A', 1), 0)
        self.assertEqual(self.cm._residue_of('A', 30), 29)

    def test_brot_phi_rotates_by_requested_angle(self):
        before = _phi(self.cm, 'A', 20)
        self.cm.apply_crot(Crot(angle='PHI', chainID='A', resid1=ResID(20), resid2=ResID(40), degrees=25.0))
        after = _phi(self.cm, 'A', 20)
        d = (after - before + 180) % 360 - 180
        self.assertAlmostEqual(abs(d), 25.0, places=3)

    def test_brot_is_rigid_on_moving_set(self):
        self.cm.apply_crot(Crot(angle='PSI', chainID='A', resid1=ResID(15), resid2=ResID(40), degrees=40.0))
        # a downstream residue moves as a rigid body: its internal C-N distance is preserved
        r = self.cm._residue_of('A', 30)
        for nm_a, nm_b in (('N', 'CA'), ('CA', 'C')):
            a0 = CoordManipulator(PSF, PDB)
            va = self.cm._res_xyz(self.cm.coords, self.cm._ridx, self.cm._names, r, nm_a)
            vb = self.cm._res_xyz(self.cm.coords, self.cm._ridx, self.cm._names, r, nm_b)
            oa = a0._res_xyz(a0.coords, a0._residue_index(), a0._names, r, nm_a)
            ob = a0._res_xyz(a0.coords, a0._ridx, a0._names, r, nm_b)
            self.assertAlmostEqual(np.linalg.norm(va - vb), np.linalg.norm(oa - ob), places=5)

    def test_brot_360_is_identity(self):
        self.cm.apply_crot(Crot(angle='PHI', chainID='A', resid1=ResID(10), resid2=ResID(30), degrees=360.0))
        self.assertTrue(np.allclose(self.cm.coords, self.orig, atol=1e-6))

    def test_chi1_moves_sidechain_not_backbone(self):
        r = self.cm._residue_of('A', 22)
        bb = self.cm.select('resid 22 and name CA')       # backbone stays
        self.cm.apply_crot(Crot(angle='CHI1', chainID='A', resid1=ResID(22), degrees=45.0))
        self.assertTrue(np.allclose(self.cm.coords[bb], self.orig[bb], atol=1e-6))
        # at least one sidechain atom moved
        sc = self.cm.select('resid 22') & ~self.cm.select('resid 22 and backbone')
        self.assertTrue(np.any(np.linalg.norm(self.cm.coords[sc] - self.orig[sc], axis=1) > 0.1))

    def test_fold_alpha_drives_backbone_to_helix(self):
        self.cm.apply_crot(Crot(angle='ALPHA', chainID='A', resid1=ResID(5),
                                resid2=ResID(15), resid3=ResID(25)))
        # interior folded residues should have near-alpha phi
        for resid in (8, 10, 12):
            self.assertAlmostEqual(_phi(self.cm, 'A', resid), -57.0, delta=1.0)

    def test_chainidmap_remaps_chain(self):
        # apply the crotation to chain A via a map that sends 'Q' -> 'A'
        cm2 = CoordManipulator(PSF, PDB)
        self.cm.apply_crot(Crot(angle='PHI', chainID='A', resid1=ResID(20), resid2=ResID(40), degrees=15.0))
        cm2.apply_crot(Crot(angle='PHI', chainID='Q', resid1=ResID(20), resid2=ResID(40), degrees=15.0),
                       chainIDmap={'Q': 'A'})
        self.assertTrue(np.allclose(self.cm.coords, cm2.coords, atol=1e-9))

    def test_scrot_chi1_moves_sidechain_beyond_cb(self):
        r = self.cm._residue_of('A', 15)
        res15 = self.cm.select('resid 15')
        # CA, CB, and the backbone atoms stay put under a chi1 rotation
        fixed = res15 & np.array([a.name in ('CA', 'CB', 'N', 'C', 'O', 'HA')
                                  for a in self.cm.atoms.data])
        self.cm.apply_scrot(1, r, 40.0)
        self.assertTrue(np.allclose(self.cm.coords[fixed], self.orig[fixed], atol=1e-6))
        # something beyond CB moved
        moved = res15 & ~fixed
        self.assertTrue(np.any(np.linalg.norm(self.cm.coords[moved] - self.orig[moved], axis=1) > 0.1))

    def test_angleijk_rejected(self):
        with self.assertRaises(CoordManipulateError):
            self.cm.apply_crot(Crot(angle='ANGLEIJK', segnamei='A', residi=ResID(1), atomi='N',
                                    segnamejk='A', residj=ResID(2), atomj='C', residk=ResID(3),
                                    atomk='O', degrees=30.0))


if __name__ == '__main__':
    unittest.main()
