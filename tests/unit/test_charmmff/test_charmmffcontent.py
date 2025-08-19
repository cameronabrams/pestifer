# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
import os
from pestifer import resources
# from pestifer.core.resourcemanager import ResourceManager
from pestifer.charmmff.charmmffcontent import CHARMMFFContent, extract_resi_pres_blocks, extract_mass_lines
from pestifer.charmmff.charmmfftop import CharmmMasses

import logging
logger = logging.getLogger(__name__)

class TestCharmmffContentCacheable(unittest.TestCase):

    def test_charmmff_content_cacheable(self):
        resource_path = os.path.dirname(resources.__file__)
        charmmff_path = os.path.join(resource_path, 'charmmff')
        C = CHARMMFFContent(charmmff_path, force_rebuild=True)
        self.assertFalse(C.from_cache)
        anotherC = CHARMMFFContent(charmmff_path)
        self.assertTrue(anotherC.from_cache)

class TestCharmmffContent(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        logging.debug("Setting up TestCharmmffContent class...")
        resource_path = os.path.dirname(resources.__file__)
        charmmff_path = os.path.join(resource_path, 'charmmff')
        cls.C = CHARMMFFContent(charmmff_path)
        # Ensure the resource manager is initialized before any tests run
        logging.debug("Done setting up TestCharmmffContent class...")

    def test_charmmff_content(self):
        self.C.clean_local_charmmff_files()
        self.assertTrue(self.C!=None)
        self.assertTrue(self.C.filenamemap!={})
        basenames=[k for k in self.C.filenamemap.keys()]
        self.assertTrue(len(basenames)>0)
        self.assertEqual(len(set(basenames)),len(basenames))  # check for duplicates
        self.assertTrue(len(self.C.streams)>0)
        self.assertEqual(self.C.streams.sort(),['prot','carb', 'na', 'lipid'].sort())
        self.assertIn('pestifer.top', self.C.custom_files)

        self.C.copy_charmmfile_local('par_all36m_prot.prm')
        self.assertTrue(os.path.exists('par_all36m_prot.prm'))
        self.C.copy_charmmfile_local('top_all36_prot.rtf')
        self.assertTrue(os.path.exists('top_all36_prot.rtf'))
        self.C.copy_charmmfile_local('toppar_water_ions.str')
        self.assertTrue(os.path.exists('toppar_water_ions.str'))
        with open('toppar_water_ions.str','r') as f:
            contents=f.read()
            self.assertTrue('! commented out by pestifer' in contents)
        self.C.copy_charmmfile_local('toppar_all36_carb_glycopeptide.str')
        self.assertTrue(os.path.exists('toppar_all36_carb_glycopeptide.str'))
        with open('toppar_all36_carb_glycopeptide.str','r') as f:
            contents=f.read()
            self.assertFalse('! commented out by pestifer' in contents)
        self.C.copy_charmmfile_local('toppar_all36_moreions.str')
        self.assertTrue(os.path.exists('toppar_all36_moreions.str'))
        self.C.clean_local_charmmff_files()
        self.assertFalse(os.path.exists('par_all36m_prot.prm'))
        self.assertFalse(os.path.exists('top_all36_prot.rtf'))
        self.assertFalse(os.path.exists('toppar_water_ions.str'))
        self.assertFalse(os.path.exists('toppar_all36_carb_glycopeptide.str'))
        self.assertFalse(os.path.exists('toppar_all36_moreions.str'))

    def test_charmmffcontent_pdbrepository_initialized(self):
        self.assertTrue(self.C!=None)
        self.assertTrue(self.C.pdbrepository!=None)
        self.assertTrue(self.C.pdbrepository.collections!=None)
        self.assertTrue(len(self.C.pdbrepository.collections)>0)
        self.assertTrue('lipid' in self.C.pdbrepository.collections)
        self.assertTrue('water_ions' in self.C.pdbrepository.collections)
        self.assertTrue('TIP3' in self.C.pdbrepository)
        self.assertTrue('PSM' in self.C.pdbrepository)

    def test_charmmffcontent_pdbrepository_checkout(self):
        c=self.C.pdbrepository.checkout('PSM')
        self.assertTrue(c!=None)
        self.assertEqual(c.get_charge(),0.0)
        params=c.get_parameters()
        self.assertIn('toppar_all36_lipid_sphingo.str',params)
        self.assertEqual(len(c.info['conformers']),10)
        self.assertAlmostEqual(c.info['conformers'][0]['head-tail-length'],26.57,places=2)
        self.assertAlmostEqual(c.info['conformers'][0]['max-internal-length'],30.54,places=2)
        c.get_pdb(0)
        self.assertTrue(os.path.exists('PSM-00.pdb'))
        os.remove('PSM-00.pdb')
        
    def test_topfile_of_patchname(self):
        """Test that the topfile of a patchname is returned correctly."""
        patchname = 'NNEU'
        topfile = self.C.get_topfile_of_patchname(patchname)
        self.assertEqual(topfile, 'top_all36_prot.rtf')
        patchname = 'TYRO'
        topfile = self.C.get_topfile_of_patchname(patchname)
        self.assertEqual(topfile, 'pestifer.top')

    def test_get_topfile_blocks(self):
        """Test that the get_topfile_blocks method returns the correct blocks."""
        topfile = 'toppar/top_all36_prot.rtf'
        content=self.C.contents_from_topfile(topfile)
        self.assertTrue(len(content) > 0)
        blocks = extract_resi_pres_blocks(content)
        self.assertTrue(isinstance(blocks, list))
        self.assertEqual(len(blocks), 48)
        self.assertTrue(all(isinstance(block, str) for block in blocks))

    def test_extract_mass_lines(self):
        """Test that the extract_mass_lines method returns the correct mass lines."""
        topfile = 'toppar/top_all36_prot.rtf'
        content=self.C.contents_from_topfile(topfile)
        self.assertTrue(len(content) > 0)
        mass_lines = extract_mass_lines(content)
        self.assertTrue(isinstance(mass_lines, list))
        self.assertTrue(all(isinstance(line, str) for line in mass_lines))
        self.assertTrue(len(mass_lines) > 0)
        self.assertEqual(len(mass_lines), 54)

    def test_resis_and_masses_from_topfile(self):
        """Test that the resis_and_masses_from_topfile method returns the correct residues and masses."""
        topfile = 'toppar/top_all36_prot.rtf'
        resis, masses = self.C.resis_and_masses_from_topfile(topfile)
        self.assertTrue(isinstance(resis, list))
        self.assertEqual(len(resis), 48)
        self.assertTrue(isinstance(masses, CharmmMasses))
        self.assertEqual(len(masses), 54)
        r=resis[0]
        self.assertTrue(r.resname == 'ALA')
