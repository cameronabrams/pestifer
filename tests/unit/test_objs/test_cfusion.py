# Author: Cameron F. Abrams, <cfa22@drexel.edu>
# with substantial contributions from Chat GPT 4o

import unittest
from pestifer.objs.cfusion import Cfusion, CfusionList

class TestCfusion(unittest.TestCase):

    def test_cfusion_creation(self):
        cfusion = Cfusion(
            sourcefile="test.pdb",
            sourceseg="A",
            resseqnum1=1,
            insertion1="",
            resseqnum2=10,
            insertion2="",
            chainID="A",
            id=1
        )
        self.assertIsInstance(cfusion, Cfusion)
        self.assertEqual(cfusion.sourcefile, "test.pdb")
        self.assertEqual(cfusion.sourceseg, "A")
        self.assertEqual(cfusion.resseqnum1, 1)
        self.assertEqual(cfusion.insertion1, "")
        self.assertEqual(cfusion.resseqnum2, 10)
        self.assertEqual(cfusion.insertion2, "")
        self.assertEqual(cfusion.chainID, "A")
        self.assertEqual(cfusion.id, 1)
        self.assertEqual(cfusion._yaml_header, 'Cfusions')
        self.assertEqual(cfusion._objcat, 'seq')
        self.assertEqual(cfusion.describe(), "Cfusion(sourcefile=test.pdb, sourceseg=A, resseqnum1=1, insertion1=, resseqnum2=10, insertion2=, chainID=A, id=1)")
        self.assertEqual(Cfusion._counter, 0)  # not created using from_input, so counter should not increment

    def test_cfusion_from_shortcode(self):
        # shortcode format: filename:C:nnn-ccc,S
        # filename -- fusion coordinate filename
        # C -- source segment
        # nnn -- N-terminal resid of fusion sequence
        # ccc -- C-terminal resid of fusion sequence
        # S -- chainID, segment in base-molecule the fusion is fused to
        shortcode = Cfusion.Adapter.from_string("test.pdb:A:1-10,A")
        cfusion = Cfusion.from_input(shortcode)
        self.assertIsInstance(cfusion, Cfusion)
        self.assertEqual(cfusion.sourcefile, "test.pdb")
        self.assertEqual(cfusion.sourceseg, "A")
        self.assertEqual(cfusion.resseqnum1, 1)
        self.assertEqual(cfusion.insertion1, "")
        self.assertEqual(cfusion.resseqnum2, 10)
        self.assertEqual(cfusion.insertion2, "")
        self.assertEqual(cfusion.chainID, "A")
        self.assertEqual(cfusion.id, 0)
        # do it again to make sure the counter increments
        cfusion2 = Cfusion.from_input(shortcode)
        self.assertEqual(cfusion2.id, 1)

    def test_cfusion_to_input_string(self):
        cfusion = Cfusion(
            sourcefile="test.pdb",
            sourceseg="A",
            resseqnum1=1,
            insertion1="",
            resseqnum2=10,
            insertion2="",
            chainID="A",
            id=1
        )
        input_string = cfusion.to_input_string()
        self.assertEqual(input_string, "test.pdb:A:1-10,A")

    def test_cfusion_list_creation(self):
        cfusion_list = CfusionList()
        self.assertIsInstance(cfusion_list, CfusionList)
        self.assertEqual(len(cfusion_list), 0)