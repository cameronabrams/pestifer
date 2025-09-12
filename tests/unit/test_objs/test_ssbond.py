import unittest

from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict

from pathlib import Path

from mmcif.api.PdbxContainers import DataContainer

from pestifer.objs.resid import ResID
from pestifer.objs.ssbond import SSBond, SSBondList

from pestifer.molecule.residue import Residue, ResidueList
from pestifer.molecule.atom import AtomList

from pestifer.psfutil.psfpatch import PSFDISUPatch

from pestifer.util.cifutil import CIFload

class TestSSBond(unittest.TestCase):

    inputs_dir = Path(__file__).parents[2] / "inputs"

    def test_ssbond_creation(self):
        ssbond = SSBond(
            chainID1="A",
            resid1=ResID(10),
            chainID2="B",
            resid2=ResID(20)
        )
        self.assertIsInstance(ssbond, SSBond)
        self.assertEqual(ssbond.chainID1, "A")
        self.assertEqual(ssbond.resid1, ResID(10))
        self.assertEqual(ssbond.residue1, None)
        self.assertEqual(ssbond.chainID2, "B")
        self.assertEqual(ssbond.resid2, ResID(20))
        self.assertEqual(ssbond.residue2, None)

    def test_ssbond_from_shortcode(self):
        ssbond = SSBond("A_10-B_20")
        self.assertIsInstance(ssbond, SSBond)
        self.assertEqual(ssbond.chainID1, "A")
        self.assertEqual(ssbond.resid1, ResID(10))
        self.assertEqual(ssbond.chainID2, "B")
        self.assertEqual(ssbond.resid2, ResID(20))

    def test_ssbond_from_pdbrecord(self):
        p = PDBParser(filepath=self.inputs_dir / '4zmj.pdb').parse()
        pr = p.parsed[SSBond._PDB_keyword][0]
        ssbond = SSBond(pr)
        self.assertIsInstance(ssbond, SSBond)

    def test_ssbond_from_psfdisupatch(self):
        patch = PSFDISUPatch(['A:10', 'B:20'])
        ssbond = SSBond(patch)
        self.assertIsInstance(ssbond, SSBond)
        self.assertEqual(ssbond.chainID1, 'A')
        self.assertEqual(ssbond.resid1, ResID(10))
        self.assertEqual(ssbond.chainID2, 'B')
        self.assertEqual(ssbond.resid2, ResID(20))

class TestSSBondList(unittest.TestCase):

    inputs_dir = Path(__file__).parents[2] / "inputs"

    def test_ssbond_list_creation(self):
        ssbond_list = SSBondList()
        self.assertIsInstance(ssbond_list, SSBondList)
        self.assertEqual(len(ssbond_list), 0)

    def test_ssbond_list_from_pdbrecords(self):
        p = PDBParser(filepath=self.inputs_dir / '4zmj.pdb').parse().parsed
        self.assertIsInstance(p, PDBRecordDict)
        ssbond_list = SSBondList.from_pdb(p)
        self.assertIsInstance(ssbond_list, SSBondList)
        self.assertGreater(len(ssbond_list), 0)

    def test_ssbond_list_from_cif(self):
        cif_data = CIFload(Path(self.inputs_dir / '4zmj.cif'))
        self.assertIsInstance(cif_data, DataContainer)
        ssbond_list = SSBondList.from_cif(cif_data)
        self.assertIsInstance(ssbond_list, SSBondList)
        self.assertGreater(len(ssbond_list), 0)

    def test_ssbond_list_set_residues(self):
        # make a mock list of 10 residues with sequential resids all in chain A
        residues = ResidueList([
            Residue(chainID="A", segname='A', resid=ResID(10), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(11), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(12), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(13), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(14), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(15), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(16), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(17), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(18), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(19), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(20), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", segname='A', resid=ResID(21), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein')
        ])
        # make a mock list of ssbonds
        ssbond_list = SSBondList([
            SSBond(chainID1="A", resid1=ResID(10), chainID2="A", resid2=ResID(20)),
            SSBond(chainID1="A", resid1=ResID(11), chainID2="A", resid2=ResID(19)),
            SSBond(chainID1="A", resid1=ResID(12), chainID2="A", resid2=ResID(18))])
        
        ssbond_list.assign_residues(residues)
        self.assertEqual(ssbond_list[0].residue1, residues[0])
        self.assertEqual(ssbond_list[0].residue2, residues[10])
        self.assertEqual(ssbond_list[1].residue1, residues[1])
        self.assertEqual(ssbond_list[1].residue2, residues[9])
        self.assertEqual(ssbond_list[2].residue1, residues[2])
        self.assertEqual(ssbond_list[2].residue2, residues[8])

