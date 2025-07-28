import unittest
from pestifer.core.objmanager import ObjManager
from pestifer.molecule.segment import Segment, SegmentList
from pestifer.molecule.chainidmanager import ChainIDManager
from pestifer.molecule.asymmetricunit import AsymmetricUnit
from pidibble.pdbparse import PDBParser

class TestAsymmetricUnit(unittest.TestCase):
    
    def test_asymmetric_unit_initialization_empty(self):
        AU = AsymmetricUnit()
        self.assertIsInstance(AU, AsymmetricUnit)
        self.assertFalse(hasattr(AU, 'segments'))

    def test_asymmetric_unit_initialization_from_pdb(self):
        p = PDBParser(filepath='fixtures/data/6pti.pdb').parse()
        AU = AsymmetricUnit(parsed=p.parsed)
        self.assertIsInstance(AU, AsymmetricUnit)
        self.assertTrue(hasattr(AU, 'segments'))
        self.assertIsInstance(AU.segments, SegmentList)