# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import unittest

from mmcif.api.PdbxContainers import DataContainer
from pathlib import Path

from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict

from pestifer.molecule.atom import AtomList
from pestifer.molecule.residue import Residue, ResidueList
from pestifer.objs.resid import ResID
from pestifer.objs.seqadv import Seqadv, SeqadvList
from pestifer.util.cifutil import CIFdict, CIFload

class TestSeqadv(unittest.TestCase):

    inputs_dir = Path(__file__).parents[2] / "inputs"

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
        p=PDBParser(filepath=self.inputs_dir / '4zmj.pdb').parse()
        pr=p.parsed[Seqadv._PDB_keyword][0]
        seqadv = Seqadv(pr)
        self.assertIsInstance(seqadv, Seqadv)

    def test_seqadv_from_cifdict(self):
        p=CIFload(self.inputs_dir / '4zmj.cif')
        obj=p.getObj(Seqadv._CIF_CategoryName)
        d=CIFdict(obj, 0)
        seqadv = Seqadv(d)
        self.assertIsInstance(seqadv, Seqadv)

class TestSeqadvList(unittest.TestCase):

    inputs_dir = Path(__file__).parents[2] / "inputs"

    def tearDown(self):
        for f in Path('.').glob('*.log'):
            f.unlink()

    def test_seqadv_list_creation(self):
        seqadv_list = SeqadvList()
        self.assertIsInstance(seqadv_list, SeqadvList)
        self.assertEqual(len(seqadv_list), 0)

    def test_seqadv_list_from_pdb(self):
        p = PDBParser(filepath=self.inputs_dir / '4zmj.pdb').parse().parsed
        self.assertIsInstance(p, PDBRecordDict)
        seqadv_list = SeqadvList.from_pdb(p)
        self.assertIsInstance(seqadv_list, SeqadvList)
        self.assertGreater(len(seqadv_list), 0)

    def test_seqadv_list_from_cif(self):
        cif_data = CIFload(self.inputs_dir / '4zmj.cif')
        self.assertIsInstance(cif_data, DataContainer)
        seqadv_list = SeqadvList.from_cif(cif_data)
        self.assertIsInstance(seqadv_list, SeqadvList)
        self.assertGreater(len(seqadv_list), 0)

    def test_seqadv_list_assign_residues(self):
        S = SeqadvList([Seqadv(idCode="1ABC", resname="ALA", chainID="A", resid=ResID(1), typekey="engineered mutation", dbRes="GLY"), 
                        Seqadv(idCode="2DEF", resname="GLY", chainID="B", resid=ResID(2), typekey="conflict", dbRes="ALA")])
        R = ResidueList([Residue(resname="ALA", chainID="A", segname='A', resid=ResID(1), segtype='protein', atoms=AtomList([]), resolved=True),
                         Residue(resname="GLY", chainID="B", segname='B',resid=ResID(2), segtype='protein', atoms=AtomList([]), resolved=True)])
        S.assign_residues(R)
        self.assertEqual(len(S), 2)
        self.assertEqual(S[0].residue, R[0])
        self.assertEqual(S[1].residue, R[1])