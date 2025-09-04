# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import unittest
from pestifer.molecule.atom import AtomList
from pestifer.molecule.residue import Residue, ResidueList
from pestifer.objs.patch import Patch, PatchList
from pestifer.objs.resid import ResID

class TestPatch(unittest.TestCase):
    def test_patch_creation(self):
        patch = Patch(
            chainID="A",
            resid=ResID(10),
            patchname="test_patch"
        )
        self.assertIsInstance(patch, Patch)
        self.assertEqual(patch.chainID, "A")
        self.assertEqual(patch.resid, ResID(10))
        self.assertTrue(patch.residue is None)
        self.assertEqual(patch.patchname, "test_patch")
        self.assertEqual(patch._yaml_header, 'patches')
        self.assertEqual(patch._objcat, 'seq')
        self.assertEqual(repr(patch), "Patch(patchname='test_patch', chainID='A', resid=ResID(resseqnum=10), use_after_regenerate=False)")

    def test_patch_from_shortcode(self):
        patch = Patch("test_patch:A:10")
        self.assertIsInstance(patch, Patch)
        self.assertEqual(patch.chainID, "A")
        self.assertEqual(patch.resid, ResID(10))
        self.assertEqual(patch.patchname, "test_patch")

    def test_patch_to_input_string(self):
        patch = Patch(
            chainID="A",
            resid=ResID(10),
            patchname="test_patch"
        )
        input_string = patch.shortcode()
        self.assertEqual(input_string, "test_patch:A:10")

    def test_patch_copy(self):
        patch = Patch(
            chainID="C",
            resid=ResID(30),
            patchname="copy_patch"
        )
        patch_copy = patch.copy()
        self.assertIsInstance(patch_copy, Patch)
        self.assertTrue(patch is not patch_copy)
        self.assertEqual(patch_copy.chainID, "C")
        self.assertEqual(patch_copy.resid, ResID(30))
        self.assertEqual(patch_copy.patchname, "copy_patch")
    
    def test_patch_list_creation(self):
        patch_list = PatchList()
        self.assertIsInstance(patch_list, PatchList)
        self.assertEqual(len(patch_list), 0)

    def test_patch_assign_residue(self):
        patch = Patch(
            chainID="A",
            resid=ResID(10),
            patchname="test_patch"
        )
        residue = Residue(
            chainID="A",
            resname="FAKEY",
            atoms=AtomList([]),
            resolved=True,
            segtype='protein',
            segname='A',
            resid=ResID(10)
        )
        patch.assign_residue(residue)
        self.assertEqual(patch.residue, residue)

    def test_patch_list_assign_residues(self):
        patch_list = PatchList([Patch(chainID="A", resid=ResID(10), patchname="patch1"), Patch(chainID="B", resid=ResID(20), patchname="patch2")])
        self.assertIsInstance(patch_list, PatchList)
        self.assertEqual(len(patch_list), 2)
        residues = ResidueList([
            Residue(chainID="A", segname='A', resname="FAKEY", atoms=AtomList([]), resolved=True, segtype='protein', resid=ResID(10)),
            Residue(chainID="B", segname='B', resname="FAKEY", atoms=AtomList([]), resolved=True, segtype='protein', resid=ResID(20))
        ])
        patch_list.assign_residues(residues)
        self.assertEqual(patch_list[0].residue, residues[0])
        self.assertEqual(patch_list[1].residue, residues[1])