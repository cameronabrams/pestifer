
import logging
import os
import unittest
from pestifer.charmmff.charmmffcontent import CHARMMFFContent
from pestifer.charmmff.charmmffresidatabase import CHARMMFFResiDatabase
from pestifer import resources

logger = logging.getLogger(__name__)

class TestCharmmffResiDatabase(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # Ensure the resource manager is initialized before any tests run
        logging.debug("Setting up TestCharmmffResiDatabase class...")
        resource_path = os.path.dirname(resources.__file__)
        charmmff_path = os.path.join(resource_path, 'charmmff')
        cls.C = CHARMMFFContent(charmmff_path)
        cls.CDB = CHARMMFFResiDatabase(cls.C,streamIDs=[])
        cls.CDB.add_stream('lipid')
        cls.CDB.add_topology('toppar_all36_moreions.str', streamIDoverride='water_ions')
        logging.debug("Done setting up TestCharmmffResiDatabase class...")

    def test_charmmff_residatabase_1(self):
        self.assertTrue(self.CDB!=None)
        self.assertEqual(len(self.CDB.residues),1496)
        self.assertEqual(len(self.CDB.patches),191)
        self.assertTrue('ALA' in self.CDB.residues)
        self.assertTrue('TIP3' in self.CDB.residues)
        self.assertTrue('FAKE' not in self.CDB.residues)
        self.assertTrue('NNEU' in self.CDB.patches)
        self.assertTrue('ASPP' in self.CDB.patches)
        self.assertTrue('GLUP' in self.CDB.patches)
        resname='ALA'
        ala=self.CDB.get_resi(resname)
        self.assertTrue(ala!=None)
        self.assertTrue(ala.metadata['streamID']=='prot')
        self.assertTrue(ala.metadata['substreamID']=='')
        self.assertEqual(ala.mass(),71.0794)
        resname='TOCL1'
        TOCL1=self.CDB.get_resi(resname)
        self.assertTrue(TOCL1!=None)
        self.assertTrue(TOCL1.metadata['streamID']=='lipid')
        self.assertTrue(TOCL1.metadata['substreamID']=='cardiolipin')
        self.assertAlmostEqual(TOCL1.mass(), 1457.02, places=2)

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

