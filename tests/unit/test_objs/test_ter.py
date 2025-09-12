import unittest
from pestifer.objs.ter import Ter, TerList
from pidibble.pdbparse import PDBParser
from pestifer.objs.resid import ResID
from pathlib import Path

class TestTer(unittest.TestCase):

    inputs_dir = Path(__file__).parents[2] / "inputs"

    def test_ter_creation(self):
        ter = Ter(
            serial=-1,
            chainID="A",
            resname="ABC",
            resid=ResID(1)
        )
        self.assertIsInstance(ter, Ter)
        self.assertEqual(ter.serial, -1)
        self.assertEqual(ter.chainID, "A")
        self.assertEqual(ter.resname, "ABC")
        self.assertEqual(ter.resid, ResID(1))
        self.assertEqual(ter._yaml_header, 'terminals')
        self.assertEqual(ter._objcat, 'seq')
        self.assertEqual(repr(ter), "Ter(serial=-1, resname='ABC', chainID='A', resid=ResID(resseqnum=1))")

    def test_ter_from_pdbrecord(self):
        p = PDBParser(filepath=self.inputs_dir / '4zmj.pdb').parse()
        pr = p.parsed[Ter._PDB_keyword][0]
        ter = Ter(pr)
        self.assertIsInstance(ter, Ter)
        self.assertEqual(ter.serial, pr.serial)
        self.assertEqual(ter.chainID, pr.residue.chainID)
        self.assertEqual(ter.resname, pr.residue.resName)
        self.assertEqual(ter.resid, ResID(resseqnum=pr.residue.seqNum, insertion=pr.residue.iCode))
