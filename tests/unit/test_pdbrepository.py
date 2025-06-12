import unittest
import os
from pestifer.pdbrepository import PDBRepository, PDBInfo, PDBCollection
import logging
logger = logging.getLogger(__name__)

class TestPDBRepository(unittest.TestCase):
    
    def test_pdbcollection_initialization(self):
        coll=PDBCollection('lipid.tgz')
        self.assertEqual(coll.streamID, 'lipid')
        self.assertEqual(coll.path, 'lipid.tgz')
        self.assertTrue('PSM' in coll)
        self.assertFalse('FAKE' in coll)
        coll2=PDBCollection('solos')
        self.assertEqual(coll2.streamID, 'solos')
        self.assertEqual(coll2.path, 'solos')
        self.assertTrue('TIP3' in coll2)
        coll3=PDBCollection('mylipid.tgz')
        self.assertEqual(coll3.streamID, 'mylipid')
        self.assertEqual(coll3.path, 'mylipid.tgz')
        self.assertTrue('C7DHPC' in coll3)

    def test_pdbrepository_initialization(self):
        pdb_repo = PDBRepository()
        self.assertIsNotNone(pdb_repo)
        pdb_repo.add_path('lipid.tgz')
        pdb_repo.add_path('solos')
        pdb_repo.add_path('mylipid.tgz')
        self.assertTrue('lipid' in pdb_repo.collections)
        self.assertTrue('solos' in pdb_repo.collections)
        self.assertTrue('mylipid' in pdb_repo.collections)
        self.assertTrue('PSM' in pdb_repo)
        self.assertTrue('TIP3' in pdb_repo)
        self.assertTrue('C7DHPC' in pdb_repo)
        self.assertFalse('FAKE' in pdb_repo)

    def test_pdbrepository_checkout_fromtar(self):
        pdb_repo = PDBRepository()
        pdb_repo.add_path('lipid.tgz')
        c = pdb_repo.checkout('PSM')
        self.assertIsNotNone(c)
        self.assertEqual(c.get_charge(), 0.0)
        params = c.get_parameters()
        self.assertIn('toppar_all36_lipid_sphingo.str', params)
        self.assertEqual(len(c.info['conformers']), 10)
        self.assertAlmostEqual(c.info['conformers'][0]['head-tail-length'], 26.779, places=3)
        self.assertAlmostEqual(c.info['conformers'][0]['max-internal-length'], 30.398, places=3)

    def test_pdbrepository_checkout_tip3(self):
        pdb_repo = PDBRepository()
        pdb_repo.add_path('solos')
        c = pdb_repo.checkout('TIP3')
        self.assertIsNotNone(c)
        self.assertEqual(c.get_charge(), 0.0)
        self.assertEqual(c.info, PDBInfo({}))
        c.get_pdb(0)
        self.assertTrue(os.path.exists('TIP3.pdb'))
        os.remove('TIP3.pdb')

    def test_pdbrepository_checkout_user(self):
        pdb_repo = PDBRepository()
        pdb_repo.add_path('mylipid.tgz')
        c = pdb_repo.checkout('C7DHPC')
        self.assertIsNotNone(c)
        c.get_pdb(1)
        self.assertTrue(os.path.exists('C7DHPC-01.pdb'))
        os.remove('C7DHPC-01.pdb')