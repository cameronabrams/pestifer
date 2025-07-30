import unittest
from pestifer.objs.patch import Patch, PatchList

class TestPatch(unittest.TestCase):
    def test_patch_creation(self):
        patch = Patch(
            chainID="A",
            resseqnum=10,
            insertion="",
            patchname="test_patch"
        )
        self.assertIsInstance(patch, Patch)
        self.assertEqual(patch.chainID, "A")
        self.assertEqual(patch.resseqnum, 10)
        self.assertEqual(patch.insertion, "")
        self.assertEqual(patch.patchname, "test_patch")
        self.assertEqual(patch._yaml_header, 'patches')
        self.assertEqual(patch._objcat, 'seq')
        self.assertEqual(repr(patch), "Patch(patchname='test_patch', chainID='A', resseqnum=10, insertion='', use_after_regenerate=False)")

    def test_patch_from_shortcode(self):
        patch = Patch("test_patch:A:10")
        self.assertIsInstance(patch, Patch)
        self.assertEqual(patch.chainID, "A")
        self.assertEqual(patch.resseqnum, 10)
        self.assertEqual(patch.insertion, "")
        self.assertEqual(patch.patchname, "test_patch")

    def test_patch_to_input_string(self):
        patch = Patch(
            chainID="A",
            resseqnum=10,
            insertion="",
            patchname="test_patch"
        )
        input_string = patch.shortcode()
        self.assertEqual(input_string, "test_patch:A:10")

    def test_patch_copy(self):
        patch = Patch(
            chainID="C",
            resseqnum=30,
            insertion="",
            patchname="copy_patch"
        )
        patch_copy = patch.copy()
        self.assertIsInstance(patch_copy, Patch)
        self.assertTrue(patch is not patch_copy)
        self.assertEqual(patch_copy.chainID, "C")
        self.assertEqual(patch_copy.resseqnum, 30)
        self.assertEqual(patch_copy.insertion, "")
        self.assertEqual(patch_copy.patchname, "copy_patch")
    
    def test_patch_list_creation(self):
        patch_list = PatchList()
        self.assertIsInstance(patch_list, PatchList)
        self.assertEqual(len(patch_list), 0)