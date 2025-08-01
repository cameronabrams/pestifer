import unittest

from pestifer.objs.ssbond import SSBond, SSBondList
from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict
from pathlib import Path
from mmcif.api.PdbxContainers import DataContainer
from pestifer.util.cifutil import CIFload
from pestifer.objs.resid import ResID
from pestifer.molecule.residue import Residue, ResidueList
from pestifer.molecule.atom import Atom, AtomList
from pestifer.objs.mutation import Mutation, MutationList
class TestSSBond(unittest.TestCase):

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
        p = PDBParser(filepath='data/4zmj.pdb').parse()
        pr = p.parsed[SSBond._PDB_keyword][0]
        ssbond = SSBond(pr)
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

    def test_ssbond_list_set_residues(self):
        # make a mock list of 10 residues with sequential resids all in chain A
        residues = ResidueList([
            Residue(chainID="A", resid=ResID(10), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(11), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(12), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(13), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(14), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(15), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(16), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(17), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(18), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(19), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(20), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein'),
            Residue(chainID="A", resid=ResID(21), resname="CYS", atoms=AtomList([]), resolved=True, segtype='protein')
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

    def test_ssbond_list_prune_muations(self):
        # make a list of 3 mock mutations
        mutations = MutationList([
            Mutation(chainID="A", origresname="CYS", resid=ResID(10), newresname="ALA", typekey="user"),
            Mutation(chainID="A", origresname="CYS", resid=ResID(11), newresname="LYS", typekey="user"),
        ])
        # make a mock list of ssbonds
        ssbond1 = SSBond(chainID1="A", resid1=ResID(10), chainID2="A", resid2=ResID(20))
        ssbond2 = SSBond(chainID1="A", resid1=ResID(11), chainID2="A", resid2=ResID(19))
        ssbond3 = SSBond(chainID1="A", resid1=ResID(12), chainID2="A", resid2=ResID(18))

        ssbond_list = SSBondList([ssbond1, ssbond2, ssbond3])

        pruned_list = ssbond_list.prune_mutations(mutations)
        self.assertEqual(len(pruned_list), 2)
        self.assertEqual(pruned_list[0].resid1, ResID(10))
        self.assertEqual(pruned_list[1].resid1, ResID(11))