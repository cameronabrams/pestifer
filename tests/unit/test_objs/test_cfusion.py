# Author: Cameron F. Abrams, <cfa22@drexel.edu>
# with substantial contributions from Chat GPT 4o

import unittest
from pestifer.objs.cfusion import Cfusion, CfusionList
from pestifer.objs.resid import ResID

class TestCfusion(unittest.TestCase):

    def test_cfusion_creation(self):
        cfusion = Cfusion(
            sourcefile="test.pdb",
            sourceseg="A",
            resid1=ResID(1),
            resid2=ResID(10),
            chainID="A"
        )
        self.assertIsInstance(cfusion, Cfusion)
        self.assertEqual(cfusion.sourcefile, "test.pdb")
        self.assertEqual(cfusion.sourceseg, "A")
        self.assertEqual(cfusion.resid1.resid, 1)
        self.assertEqual(cfusion.resid2.resid, 10)
        self.assertEqual(cfusion.chainID, "A")
        self.assertEqual(cfusion.obj_id, 0)
        self.assertEqual(cfusion._yaml_header, 'Cfusions')
        self.assertEqual(cfusion._objcat, 'seq')
        self.assertEqual(repr(cfusion), "Cfusion(sourcefile='test.pdb', sourceseg='A', resid1=ResID(resseqnum=1), resid2=ResID(resseqnum=10), chainID='A', obj_id=0)")
        self.assertEqual(Cfusion._counter, 0)  # not created using from_input, so counter should not increment

    def test_cfusion_from_shortcode(self):
        # shortcode format: filename:C:nnn-ccc,S
        # filename -- fusion coordinate filename
        # C -- source segment
        # nnn -- N-terminal resid of fusion sequence
        # ccc -- C-terminal resid of fusion sequence
        # S -- chainID, segment in base-molecule the fusion is fused to
        # shortcode = Cfusion.Adapter.from_string("test.pdb:A:1-10,A")
        # cfusion = Cfusion.from_input(shortcode)
        cfusion = Cfusion("test.pdb:A:1-10,A")
        self.assertIsInstance(cfusion, Cfusion)
        self.assertEqual(cfusion.sourcefile, "test.pdb")
        self.assertEqual(cfusion.sourceseg, "A")
        self.assertEqual(cfusion.resid1.resid, 1)
        self.assertEqual(cfusion.resid2.resid, 10)
        self.assertEqual(cfusion.chainID, "A")
        self.assertEqual(cfusion.obj_id, 0)
        # do it again to make sure the counter increments
        cfusion2 = Cfusion("test.pdb:A:1-10,A")
        self.assertEqual(cfusion2.obj_id, 1)

    def test_cfusion_copy(self):
        cfusion1 = Cfusion(
            sourcefile="test.pdb",
            sourceseg="A",
            resid1=ResID(1),
            resid2=ResID(10),
            chainID="A"
        )
        cfusion2 = cfusion1.copy()
        self.assertTrue(cfusion1 is not cfusion2)
        self.assertEqual(cfusion1, cfusion2)
        self.assertIsInstance(cfusion2, Cfusion)
        self.assertEqual(cfusion2.sourcefile, "test.pdb")
        self.assertEqual(cfusion2.sourceseg, "A")
        self.assertEqual(cfusion2.resid1.resid, 1)
        self.assertEqual(cfusion2.resid2.resid, 10)
        self.assertEqual(cfusion2.chainID, "A")
        self.assertEqual(cfusion2.obj_id, 0)

    def test_cfusion_to_input_string(self):
        cfusion = Cfusion(
            sourcefile="test.pdb",
            sourceseg="A",
            resid1=ResID(1),
            resid2=ResID(10),
            chainID="A",
            obj_id=1
        )
        input_string = cfusion.shortcode()
        self.assertEqual(input_string, "test.pdb:A:1-10,A")

    def test_cfusion_list_creation(self):
        cfusion_list = CfusionList()
        self.assertIsInstance(cfusion_list, CfusionList)
        self.assertEqual(len(cfusion_list), 0)