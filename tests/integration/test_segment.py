import unittest
from pestifer.molecule import Molecule
from pestifer.config import Config
from pestifer.chainids import ChainIDManager
class TestSegment(unittest.TestCase):

    def test_segment(self):
        self.config=ConfigSetup('')
        chainIDmanager=ChainIDManager()
        m=Molecule(source='4zmj').activate_biological_assembly(1,chainIDmanager)
        au=m.asymmetric_unit
        protein_segments=[x for x in au.Segments if x.segtype=='PROTEIN']
        self.assertTrue(len(protein_segments),2)
        expected_set_of_chainIDs={'G','B','A','D','C'}
        self.assertEqual(set(au.chainIDs),expected_set_of_chainIDs)


 