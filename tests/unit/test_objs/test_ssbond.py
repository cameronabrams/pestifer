import unittest

from pestifer.objs.ssbond import SSBond, SSBondList
from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict
from pathlib import Path
from mmcif.api.PdbxContainers import DataContainer
from pestifer.util.cifutil import CIFload
class TestSSBond(unittest.TestCase):

    def test_ssbond_creation(self):
        ssbond = SSBond(
            chainID1="A",
            resseqnum1=10,
            insertion1="",
            chainID2="B",
            resseqnum2=20,
            insertion2=""
        )
        self.assertIsInstance(ssbond, SSBond)
        self.assertEqual(ssbond.chainID1, "A")
        self.assertEqual(ssbond.resseqnum1, 10)
        self.assertEqual(ssbond.insertion1, "")
        self.assertEqual(ssbond.chainID2, "B")
        self.assertEqual(ssbond.resseqnum2, 20)
        self.assertEqual(ssbond.insertion2, "")

    def test_ssbond_from_shortcode(self):
        ssbond = SSBond.new("A_10-B_20")
        self.assertIsInstance(ssbond, SSBond)
        self.assertEqual(ssbond.chainID1, "A")
        self.assertEqual(ssbond.resseqnum1, 10)
        self.assertEqual(ssbond.insertion1, "")
        self.assertEqual(ssbond.chainID2, "B")
        self.assertEqual(ssbond.resseqnum2, 20)
        self.assertEqual(ssbond.insertion2, "")

    def test_ssbond_from_pdbrecord(self):
        p = PDBParser(filepath='data/4zmj.pdb').parse()
        pr = p.parsed[SSBond._PDB_keyword][0]
        ssbond = SSBond.new(pr)
        self.assertIsInstance(ssbond, SSBond)

class TestSSBondList(unittest.TestCase):

    def test_ssbond_list_creation(self):
        ssbond_list = SSBondList()
        self.assertIsInstance(ssbond_list, SSBondList)
        self.assertEqual(len(ssbond_list), 0)

    def test_ssbond_list_from_pdbrecords(self):
        p = PDBParser(filepath='data/4zmj.pdb').parse().parsed
        self.assertIsInstance(p, PDBRecordDict)
        ssbond_list = SSBondList.from_pdb(p)
        self.assertIsInstance(ssbond_list, SSBondList)
        self.assertGreater(len(ssbond_list), 0)

    def test_ssbond_list_from_cif(self):
        cif_data = CIFload(Path('data/4zmj.cif'))
        self.assertIsInstance(cif_data, DataContainer)
        ssbond_list = SSBondList.from_cif(cif_data)
        self.assertIsInstance(ssbond_list, SSBondList)
        self.assertGreater(len(ssbond_list), 0)