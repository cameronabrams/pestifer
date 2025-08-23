import unittest
import os
from pestifer import resources
from pestifer.charmmff.pdbrepository import PDBInput, PDBRepository, PDBCollection
import logging
logger = logging.getLogger(__name__)

class TestPDBCollection(unittest.TestCase):

    def test_pdbcollection_initialization(self):
        coll = PDBCollection.build_from_resources('lipid.tgz')
        self.assertEqual(coll.streamID, 'lipid')
        self.assertEqual(coll.path_or_tarball, 'lipid.tgz')
        self.assertTrue('PSM' in coll)
        self.assertFalse('FAKE' in coll)
        coll2 = PDBCollection.build_from_resources('solos')
        self.assertEqual(coll2.streamID, 'solos')
        self.assertEqual(coll2.path_or_tarball, 'solos')
        self.assertTrue('TIP3' in coll2)
        coll3 = PDBCollection.build_from_resources('mylipid.tgz')
        self.assertEqual(coll3.streamID, 'mylipid')
        self.assertEqual(coll3.path_or_tarball, 'mylipid.tgz')
        self.assertTrue('C7DHPC' in coll3)

class TestPDBRepository(unittest.TestCase):
    
    def setUp(self):
        charmmffpath = os.path.join(os.path.dirname(resources.__file__),'charmmff')
        self.repopath = os.path.join(charmmffpath, 'pdbrepository')
        self.pdbrepo = PDBRepository(self.repopath)

    def test_pdbrepository_force_rebuild(self):
        self.fresh_repo = PDBRepository(self.repopath, force_rebuild=True)
        self.assertFalse(self.fresh_repo.from_cache)

    def test_pdbrepository_initialization(self):
        self.pdbrepo.add_resource('solos')
        self.pdbrepo.add_resource('mylipid.tgz')
        self.assertTrue('lipid' in self.pdbrepo.collections)
        self.assertTrue('solos' in self.pdbrepo.collections)
        self.assertTrue('mylipid' in self.pdbrepo.collections)
        self.assertTrue('PSM' in self.pdbrepo)
        self.assertTrue('TIP3' in self.pdbrepo)
        self.assertTrue('C7DHPC' in self.pdbrepo)
        self.assertFalse('FAKE' in self.pdbrepo)

    def test_pdbrepository_restricted(self):
        resnames = ['PSM', 'CHL1', 'DOPC']
        self.fresh_repo = PDBRepository(self.repopath, resnames=resnames)
        self.assertTrue('lipid' in self.fresh_repo.collections)
        self.assertTrue('PSM' in self.fresh_repo)
        self.assertTrue('CHL1' in self.fresh_repo)
        self.assertTrue('DOPC' in self.fresh_repo)
        self.assertFalse('TIP3' in self.fresh_repo)

    def test_pdbrepository_checkout_fromtar(self):
        c = self.pdbrepo.checkout('PSM')
        self.assertIsInstance(c, PDBInput)
        self.assertIsNotNone(c)
        self.assertEqual(c.get_charge(), 0.0)
        params = c.get_parameters()
        self.assertIn('toppar_all36_lipid_sphingo.str', params)
        self.assertEqual(len(c.info['conformers']), 10)
        self.assertAlmostEqual(c.info['conformers'][0]['head-tail-length'], 26.57, places=3)
        self.assertAlmostEqual(c.info['conformers'][0]['max-internal-length'], 30.54, places=3)

    def test_pdbrepository_checkout_tip3(self):
        c = self.pdbrepo.checkout('TIP3')
        self.assertIsNotNone(c)
        self.assertEqual(c.get_charge(), 0.0)
        self.assertEqual(c.info, {})
        c.get_pdb(0)
        self.assertTrue(os.path.exists('TIP3.pdb'))
        os.remove('TIP3.pdb')

    def test_pdbrepository_checkout_user(self):
        self.pdbrepo.add_resource('mylipid.tgz')
        c = self.pdbrepo.checkout('C7DHPC')
        self.assertIsNotNone(c)
        c.get_pdb(1)
        self.assertTrue(os.path.exists('C7DHPC-01.pdb'))
        os.remove('C7DHPC-01.pdb')