import unittest
import os
from pathlib import Path
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

    def test_pdbcollection_box_and_molecule_kinds(self):
        # a collection with a kind:box entry (solvent box) alongside a solo molecule entry;
        # the loader must distinguish them and expose the box's psf/pdb/edge/key-atom
        import tempfile
        import yaml
        d = tempfile.mkdtemp()
        coll_dir = os.path.join(d, 'testsolvent')
        os.makedirs(os.path.join(coll_dir, 'GROL'))
        with open(os.path.join(coll_dir, 'GROL', 'GROL-box.pdb'), 'w') as f:
            f.write('ATOM      1  C1  GROL X   1       0.0     0.0     0.0  1.00  0.00      X    C\nEND\n')
        with open(os.path.join(coll_dir, 'GROL', 'GROL-box.psf'), 'w') as f:
            f.write('PSF\n')
        with open(os.path.join(coll_dir, 'GROL', 'info.yaml'), 'w') as f:
            yaml.safe_dump({'kind': 'box', 'resname': 'GROL', 'box_edge': 24.13,
                            'key_atom': 'C1', 'nmol': 1, 'density': 1.26,
                            'psf': 'GROL-box.psf', 'pdb': 'GROL-box.pdb'}, f)
        with open(os.path.join(coll_dir, 'SOD.pdb'), 'w') as f:  # solo -> kind defaults to molecule
            f.write('ATOM      1 SOD  SOD X   1       0.0     0.0     0.0  1.00  0.00      X   NA\nEND\n')
        coll = PDBCollection.build_from_resources(coll_dir)

        box = coll.contents['GROL']
        self.assertTrue(box.is_box())
        self.assertEqual(box.get_box_edge(), 24.13)
        self.assertEqual(box.get_key_atom(), 'C1')
        cwd = os.getcwd()
        os.chdir(d)
        try:
            self.assertEqual(box.get_box_psf(), 'GROL-box.psf')
            self.assertEqual(box.get_box_pdb(), 'GROL-box.pdb')
            self.assertTrue(os.path.isfile('GROL-box.psf') and os.path.isfile('GROL-box.pdb'))
        finally:
            os.chdir(cwd)

        mol = coll.contents['SOD']   # a molecule entry is not a box
        self.assertFalse(mol.is_box())
        self.assertIsNone(mol.get_box_edge())
        self.assertIsNone(mol.get_box_psf())

class TestPDBRepository(unittest.TestCase):
    
    def setUp(self):
        charmmff_root = Path(resources.__file__).parent / 'charmmff'
        # pick the most recent charmmff *version* dir that actually holds a pdbrepository --
        # not merely the newest subdir (e.g. 'custom/', which has no pdbrepository, can be
        # touched more recently and would otherwise be chosen, breaking these tests).
        version_dirs = sorted((p for p in charmmff_root.iterdir() if (p / 'pdbrepository').is_dir()),
                              key=lambda p: p.stat().st_mtime)
        charmmff_path = version_dirs[-1]
        self.repopath = str(charmmff_path / 'pdbrepository')
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