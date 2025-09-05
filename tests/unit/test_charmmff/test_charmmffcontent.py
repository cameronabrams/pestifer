# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
import os
from pestifer import resources
from pestifer.charmmff.charmmffcontent import CHARMMFFContent

import logging
logger = logging.getLogger(__name__)
from pathlib import Path

class TestCharmmffContent(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        logging.debug("Setting up TestCharmmffContent class...")
        resource_path = Path(resources.__file__).parent
        charmmff_path = resource_path / 'charmmff'
        cls.C = CHARMMFFContent(charmmff_path)
        logging.debug("Done setting up TestCharmmffContent class...")

    def test_charmmffcontent_full_provisioning(self):
        self.C.deprovision()
        self.assertEqual(len(self.C.residues), 0)
        self.assertEqual(len(self.C.patches), 0)
        self.C.provision()
        self.assertEqual(len(self.C.residues), 2474)
        self.assertEqual(len(self.C.patches), 790)
        self.assertEqual(len(self.C.pdbrepository.collections), 2)
        self.assertEqual(len(self.C.pdbrepository.collections['lipid'].info), 130)
        self.assertEqual(len(self.C.pdbrepository.collections['water_ions'].info), 12)

    def test_charmmffcontent_restricted_provisioning(self):
        self.C.deprovision()
        self.assertEqual(len(self.C.residues), 0)
        self.assertEqual(len(self.C.patches), 0)
        self.C.provision(resnames=['ALA', 'TIP3', 'PSM', 'CHL1', 'DOPC', 'NNEU'])
        self.assertEqual(len(self.C.residues), 5)
        self.assertEqual(len(self.C.patches), 1)
        self.assertTrue('ALA' in self.C.residues)
        self.assertFalse('LYS' in self.C.residues)
        self.assertTrue('TIP3' in self.C.residues)
        self.assertTrue('TIP3' in self.C.pdbrepository)
        self.assertTrue('PSM' in self.C.residues)
        self.assertTrue('PSM' in self.C.pdbrepository)
        self.assertTrue('CHL1' in self.C.residues)
        self.assertTrue('DOPC' in self.C.residues)
        self.assertTrue('NNEU' in self.C.patches)

    def test_charmmffcontent_initialized(self):
        self.C.clean_local_charmmff_files()
        self.assertTrue(self.C != None)
        self.assertTrue(self.C.filenamemap != {})
        basenames = [k for k in self.C.filenamemap['top'].keys()]
        basenames.extend([k for k in self.C.filenamemap['toppar'].keys()])
        basenames.extend([k for k in self.C.filenamemap['par'].keys()])
        self.assertEqual(len(basenames), 54)
        self.assertEqual(len(set(basenames)), len(basenames))  # check for duplicates
        self.assertTrue(len(self.C.streams) > 0)
        self.assertEqual(self.C.streams.sort(), ['prot', 'carb', 'na', 'lipid'].sort())
        # self.assertIn('pestifer.top', self.C.custom_files)
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
        self.C.provision()
        self.assertTrue(self.C.pdbrepository!=None)
        self.assertTrue(self.C.pdbrepository.collections!=None)
        self.assertTrue(len(self.C.pdbrepository.collections)>0)
        self.assertTrue('lipid' in self.C.pdbrepository.collections)
        self.assertTrue('water_ions' in self.C.pdbrepository.collections)
        self.assertTrue('TIP3' in self.C.pdbrepository)
        self.assertTrue('PSM' in self.C.pdbrepository)

    def test_charmmffcontent_pdbrepository_checkout(self):
        self.C.provision()
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

    def test_charmffcontent_get_topfile_of_patchname(self):
        """Test that the topfile of a patchname is returned correctly."""
        patchname = 'NNEU'
        topfile = self.C.resi_to_topfile_map[patchname]
        self.assertEqual(topfile, 'top_all36_prot.rtf')
        patchname = 'TYRO'
        topfile = self.C.resi_to_topfile_map[patchname]
        self.assertEqual(topfile, 'pestifer.top')

    def test_charmmffcontent_get_charge(self):
        self.C.provision()
        resname='POT'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi!=None)
        self.assertTrue(resi.metadata['streamID']=='water_ions')
        self.assertEqual(resi.charge,1.0)
        resname='CLA'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi!=None)
        self.assertTrue(resi.metadata['streamID']=='water_ions')
        self.assertEqual(resi.charge,-1.0)
        resname='PO4'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi!=None)
        self.assertTrue(resi.metadata['streamID']=='water_ions')
        self.assertEqual(resi.charge,-3.0)

    def test_charmmffcontent_get_resi(self):
        self.C.provision()
        self.assertEqual(len(self.C.residues),2474)
        self.assertEqual(len(self.C.patches),790)
        self.assertTrue('ALA' in self.C.residues)
        self.assertTrue('TIP3' in self.C.residues)
        self.assertTrue('FAKE' not in self.C.residues)
        self.assertTrue('NNEU' in self.C.patches)
        self.assertTrue('ASPP' in self.C.patches)
        self.assertTrue('GLUP' in self.C.patches)
        resname='ALA'
        ala=self.C.get_resi(resname)
        self.assertTrue(ala!=None)
        self.assertTrue(ala.metadata['streamID']=='prot')
        self.assertTrue(ala.metadata['substreamID']=='')
        self.assertAlmostEqual(ala.mass,71.0794,places=3)
        resname='TOCL1'
        TOCL1=self.C.get_resi(resname)
        self.assertTrue(TOCL1!=None)
        self.assertTrue(TOCL1.metadata['streamID']=='lipid')
        self.assertTrue(TOCL1.metadata['substreamID']=='cardiolipin')
        self.assertAlmostEqual(TOCL1.mass, 1457.02, places=2)

    def test_charmffcontent_detect_structure_CHL1(self):
        self.C.provision()
        resname='CHL1'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='cholesterol')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertTrue(heads==['C3'])
        self.assertTrue(tails==['C27'])

    def test_charmffcontent_detect_structure_CHM1(self):
        self.C.provision()
        resname='CHM1'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='cholesterol')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertTrue(heads==['C1'])
        self.assertTrue(tails==['C6'])

    def test_charmffcontent_detect_structure_DPPC(self):
        self.C.provision()
        resname='DPPC'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        shortest_paths=resi.annotation['shortest_paths']
        dist1=shortest_paths[heads[0]][tails[0]]
        dist2=shortest_paths[heads[0]][tails[1]]
        # logger.debug(f'dist1 {dist1} dist2 {dist2}')
        self.assertEqual(dist1,25)
        self.assertEqual(dist2,26)

    def test_charmffcontent_detect_structure_SDS(self):
        self.C.provision()
        resname='SDS'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='detergent')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertTrue(heads==['S'])
        self.assertTrue(tails==['C12'])

    def test_charmffcontent_detect_structure_TOCL1(self):
        self.C.provision()
        resname='TOCL1'
        resi=self.C.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='cardiolipin')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertEqual(len(tails),4)
        self.assertTrue(heads==['C2'])
        self.assertTrue(tails==['CA18', 'CB18', 'CC18', 'CD18'])