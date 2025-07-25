import unittest
from pestifer.objs.ter import Ter, TerList
from pidibble.pdbparse import PDBParser

class TestTer(unittest.TestCase):
    def test_ter_creation(self):
        ter = Ter(
            serial=-1,
            chainID="A",
            resname="ABC",
            resseqnum=1,
            insertion=""
        )
        self.assertIsInstance(ter, Ter)
        self.assertEqual(ter.serial, -1)
        self.assertEqual(ter.chainID, "A")
        self.assertEqual(ter.resname, "ABC")
        self.assertEqual(ter.resseqnum, 1)
        self.assertEqual(ter.insertion, "")
        self.assertEqual(ter._yaml_header, 'terminals')
        self.assertEqual(ter._objcat, 'seq')
        self.assertEqual(ter.describe(), "Ter(serial=-1, chainID=A, resname=ABC, resseqnum=1, insertion=)")
    
    def test_ter_from_pdbrecord(self):
        p = PDBParser(filepath='data/4zmj.pdb').parse()
        pr = p.parsed[Ter._PDB_keyword][0]
        ter = Ter.new(pr)
        self.assertIsInstance(ter, Ter)
        self.assertEqual(ter.serial, pr.serial)
        self.assertEqual(ter.chainID, pr.residue.chainID)
        self.assertEqual(ter.resname, pr.residue.resName)
        self.assertEqual(ter.resseqnum, pr.residue.seqNum)
        self.assertEqual(ter.insertion, pr.residue.iCode)
