import unittest
import os
from pestifer.pdbrepository import PDBRepository
import logging
logger = logging.getLogger(__name__)

class TestPDBRepository(unittest.TestCase):
    
    def test_pdbrepository_initialization(self):
        pdb_repo = PDBRepository(basepath='.')
        self.assertIsNotNone(pdb_repo)
        self.assertTrue('base' in pdb_repo.collections)
        coll=pdb_repo.collections['base']
        self.assertTrue('streams' in coll)
        self.assertEqual(coll['path'], '.')
        self.assertTrue('lipid' in coll['streams'])
        self.assertTrue('solos' in coll['streams'])
        self.assertTrue('resnames' in coll['streams']['lipid'])
        self.assertTrue('PSM' in coll['streams']['lipid']['resnames'])

    def test_pdbrepository_checkout_fromtar(self):
        pdb_repo = PDBRepository(basepath='.')
        c = pdb_repo.checkout('PSM')
        self.assertIsNotNone(c)
        self.assertEqual(c.get_charge(), 0.0)
        params = c.get_parameters()
        self.assertIn('toppar_all36_lipid_sphingo.str', params)
        self.assertEqual(len(c.info['conformers']), 10)
        self.assertAlmostEqual(c.info['conformers'][0]['head-tail-length'], 26.779, places=3)
        self.assertAlmostEqual(c.info['conformers'][0]['max-internal-length'], 30.398, places=3)

    def test_pdbrepository_checkout_tip3(self):
        pdb_repo = PDBRepository(basepath='.')
        c = pdb_repo.checkout('TIP3')
        self.assertIsNotNone(c)
        self.assertEqual(c.get_charge(), 0.0)
        self.assertEqual(c.info, {})
        c.get_pdb(0)
        self.assertTrue(os.path.exists('TIP3.pdb'))
        os.remove('TIP3.pdb')

    def test_pdbrepository_checkout_user(self):
        pdb_repo = PDBRepository(basepath='.')
        pdb_repo.registercollection('users', 'user')
        c = pdb_repo.checkout('C7DHPC')
        self.assertIsNotNone(c)
        c.get_pdb(1)
        self.assertTrue(os.path.exists('C7DHPC-01.pdb'))
        os.remove('C7DHPC-01.pdb')