import unittest
from pestifer.objs.seqadv import Seqadv, SeqadvList
from pidibble.pdbparse import PDBParser
from pestifer.util.cifutil import CIFdict, CIFload
from pathlib import Path

class TestSeqadv(unittest.TestCase):

    def test_seqadv_creation(self):
        seqadv = Seqadv(
            idCode="1ABC",
            resname="ALA",
            chainID="A",
            resseqnum=1,
            insertion="",
            typekey="conflict"
        )
        self.assertIsInstance(seqadv, Seqadv)

    def test_seqadv_from_pdbrecord(self):
        p=PDBParser(filepath='data/4zmj.pdb').parse()
        pr=p.parsed[Seqadv._PDB_keyword][0]
        seqadv = Seqadv.new(pr)
        self.assertIsInstance(seqadv, Seqadv)

    def test_seqadv_from_cifdict(self):
        p=CIFload(Path('data/4zmj.cif'))
        obj=p.getObj(Seqadv._mmCIF_name)
        d=CIFdict(obj, 0)
        seqadv = Seqadv.new(d)
        self.assertIsInstance(seqadv, Seqadv)

    def test_seqadv_list_creation(self):
        seqadv_list = SeqadvList()
        self.assertIsInstance(seqadv_list, SeqadvList)
        self.assertEqual(len(seqadv_list), 0)