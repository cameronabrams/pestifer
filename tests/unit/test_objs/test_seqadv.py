import unittest
from pestifer.objs.seqadv import Seqadv, SeqadvList
from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict
from pestifer.util.cifutil import CIFdict, CIFload
from pathlib import Path
from mmcif.api.PdbxContainers import DataContainer
from pestifer.objs.resid import ResID
class TestSeqadv(unittest.TestCase):

    def test_seqadv_creation(self):
        seqadv = Seqadv(
            idCode="1ABC",
            resname="ALA",
            chainID="A",
            resid=ResID(1),
            typekey="engineered mutation",
            dbRes="GLY"
        )
        self.assertIsInstance(seqadv, Seqadv)
        self.assertEqual(repr(seqadv), "Seqadv(idCode='1ABC', resname='ALA', chainID='A', resid=ResID(resseqnum=1), dbRes='GLY', typekey='engineered mutation')")

    def test_seqadv_from_pdbrecord(self):
        p=PDBParser(filepath='data/4zmj.pdb').parse()
        pr=p.parsed[Seqadv._PDB_keyword][0]
        seqadv = Seqadv(pr)
        self.assertIsInstance(seqadv, Seqadv)

    def test_seqadv_from_cifdict(self):
        p=CIFload(Path('data/4zmj.cif'))
        obj=p.getObj(Seqadv._CIF_CategoryName)
        d=CIFdict(obj, 0)
        seqadv = Seqadv(d)
        self.assertIsInstance(seqadv, Seqadv)

class TestSeqadvList(unittest.TestCase):

    def test_seqadv_list_creation(self):
        seqadv_list = SeqadvList()
        self.assertIsInstance(seqadv_list, SeqadvList)
        self.assertEqual(len(seqadv_list), 0)

    def test_seqadv_list_from_pdb(self):
        p = PDBParser(filepath='data/4zmj.pdb').parse().parsed
        self.assertIsInstance(p, PDBRecordDict)
        seqadv_list = SeqadvList.from_pdb(p)
        self.assertIsInstance(seqadv_list, SeqadvList)
        self.assertGreater(len(seqadv_list), 0)

    def test_seqadv_list_from_cif(self):
        cif_data = CIFload(Path('data/4zmj.cif'))
        self.assertIsInstance(cif_data, DataContainer)
        seqadv_list = SeqadvList.from_cif(cif_data)
        self.assertIsInstance(seqadv_list, SeqadvList)
        self.assertGreater(len(seqadv_list), 0)