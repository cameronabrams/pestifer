# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
import os
from pestifer.resourcemanager import ResourceManager
from pestifer.charmmffcontent import CHARMMFFResiDatabase

import logging
logger = logging.getLogger(__name__)

class TestCharmmffContent(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Ensure the resource manager is initialized before any tests run
        logging.debug("Setting up TestCharmmffContent class...")
        cls.RM = ResourceManager()
        cls.C = cls.RM.charmmff_content
        cls.CDB = CHARMMFFResiDatabase(cls.C,streamIDs=[])
        cls.CDB.add_stream('lipid')
        cls.CDB.add_topology('toppar_all36_moreions.str', streamIDoverride='water_ions')
        logging.debug("Done setting up TestCharmmffContent class...")

    @classmethod
    def tearDownClass(cls):
        del cls.RM.charmmff_content

    def test_charmmff_content(self):
        self.RM.charmmff_content.clean_local_charmmff_files()
        self.assertTrue(self.RM.charmmff_content!=None)
        self.assertTrue(self.RM.charmmff_content.filenamemap!={})
        basenames=[k for k in self.RM.charmmff_content.filenamemap.keys()]
        self.assertTrue(len(basenames)>0)
        self.assertEqual(len(set(basenames)),len(basenames))  # check for duplicates
        self.assertTrue(len(self.RM.charmmff_content.streams)>0)
        self.assertEqual(self.RM.charmmff_content.streams,['prot','carb', 'na', 'lipid'])
        self.assertTrue(self.RM.charmmff_content.custom_files!=None)

        self.RM.charmmff_content.copy_charmmfile_local('par_all36m_prot.prm')
        self.assertTrue(os.path.exists('par_all36m_prot.prm'))
        self.RM.charmmff_content.copy_charmmfile_local('top_all36_prot.rtf')
        self.assertTrue(os.path.exists('top_all36_prot.rtf'))
        self.RM.charmmff_content.copy_charmmfile_local('toppar_water_ions.str')
        self.assertTrue(os.path.exists('toppar_water_ions.str'))
        with open('toppar_water_ions.str','r') as f:
            contents=f.read()
            self.assertTrue('! commented out by pestifer' in contents)
        self.RM.charmmff_content.copy_charmmfile_local('toppar_all36_carb_glycopeptide.str')
        self.assertTrue(os.path.exists('toppar_all36_carb_glycopeptide.str'))
        with open('toppar_all36_carb_glycopeptide.str','r') as f:
            contents=f.read()
            self.assertFalse('! commented out by pestifer' in contents)
        self.RM.charmmff_content.copy_charmmfile_local('toppar_all36_moreions.str')
        self.assertTrue(os.path.exists('toppar_all36_moreions.str'))
        self.RM.charmmff_content.clean_local_charmmff_files()
        self.assertFalse(os.path.exists('par_all36m_prot.prm'))
        self.assertFalse(os.path.exists('top_all36_prot.rtf'))
        self.assertFalse(os.path.exists('toppar_water_ions.str'))
        self.assertFalse(os.path.exists('toppar_all36_carb_glycopeptide.str'))
        self.assertFalse(os.path.exists('toppar_all36_moreions.str'))

    def test_charmmff_residatabase_1(self):
        self.assertTrue(self.CDB!=None)
        self.assertEqual(len(self.CDB.residues),1496)
        self.assertTrue('ALA' in self.CDB.residues)
        self.assertTrue('TIP3' in self.CDB.residues)
        self.assertTrue('FAKE' not in self.CDB.residues)
        resname='ALA'
        ala=self.CDB.get_resi(resname)
        self.assertTrue(ala!=None)
        self.assertTrue(ala.metadata['streamID']=='prot')
        self.assertTrue(ala.metadata['substreamID']=='')
        self.assertEqual(ala.mass(),71.0794)

    def test_detect_charge(self):
        resname='POT'
        resi=self.CDB.get_resi(resname)
        self.assertTrue(resi!=None)
        self.assertTrue(resi.metadata['streamID']=='water_ions')
        self.assertEqual(resi.charge,1.0)
        resname='CLA'
        resi=self.CDB.get_resi(resname)
        self.assertTrue(resi!=None)
        self.assertTrue(resi.metadata['streamID']=='water_ions')
        self.assertEqual(resi.charge,-1.0)
        resname='PO4'
        resi=self.CDB.get_resi(resname)
        self.assertTrue(resi!=None)
        self.assertTrue(resi.metadata['streamID']=='water_ions')
        self.assertEqual(resi.charge,-3.0)

    def test_detect_structure_CHL1(self):
        resname='CHL1'
        resi=self.CDB.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='cholesterol')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertTrue(heads==['C3'])
        self.assertTrue(tails==['C27'])

    def test_detect_structure_CHM1(self):
        resname='CHM1'
        resi=self.CDB.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='cholesterol')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertTrue(heads==['C1'])
        self.assertTrue(tails==['C6'])

    def test_detect_structure_DPPC(self):
        resname='DPPC'
        resi=self.CDB.get_resi(resname)
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

    def test_detect_structure_SDS(self):
        resname='SDS'
        resi=self.CDB.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='detergent')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertTrue(heads==['S'])
        self.assertTrue(tails==['C12'])

    def test_detect_structure_TOCL1(self):
        resname='TOCL1'
        resi=self.CDB.get_resi(resname)
        self.assertTrue(resi.metadata['substreamID']=='cardiolipin')
        resi.lipid_annotate()
        heads=resi.annotation['heads']
        tails=resi.annotation['tails']
        self.assertEqual(len(tails),4)
        self.assertTrue(heads==['C2'])
        self.assertTrue(tails==['CA18', 'CB18', 'CC18', 'CD18'])

    def test_charmmffcontent_pdbrepository_initialized(self):
        self.assertTrue(self.RM!=None)
        self.assertTrue(self.RM.charmmff_content!=None)
        self.assertTrue(self.RM.pdbrepository!=None)
        self.assertTrue(self.RM.pdbrepository.collections!=None)
        self.assertTrue(len(self.RM.pdbrepository.collections)>0)
        self.assertTrue('lipid' in self.RM.pdbrepository.collections)
        self.assertTrue('water_ions' in self.RM.pdbrepository.collections)
        self.assertTrue('TIP3' in self.RM.pdbrepository)
        self.assertTrue('PSM' in self.RM.pdbrepository)

    def test_charmmffcontent_pdbrepository_checkout(self):
        c=self.RM.pdbrepository.checkout('PSM')
        self.assertTrue(c!=None)
        self.assertEqual(c.get_charge(),0.0)
        params=c.get_parameters()
        self.assertIn('toppar_all36_lipid_sphingo.str',params)
        self.assertEqual(len(c.info['conformers']),10)
        self.assertAlmostEqual(c.info['conformers'][0]['head-tail-length'],27.139,places=2)
        self.assertAlmostEqual(c.info['conformers'][0]['max-internal-length'],30.72,places=2)
        c.get_pdb(0)
        self.assertTrue(os.path.exists('PSM-00.pdb'))
        os.remove('PSM-00.pdb')
        