import unittest
from pestifer.objs.graft import Graft, GraftList

class TestGraft(unittest.TestCase):
    def test_graft_creation(self):
        graft = Graft(
            orig_chainID="A",
            orig_resseqnum1=520,
            orig_insertion1="",
            source_pdbid="4ijk",
            source_chainID="A",
            source_resseqnum1=1,
            source_insertion1="",
        )
        self.assertIsInstance(graft, Graft)
        self.assertEqual(graft.orig_chainID, "A")
        self.assertEqual(graft.orig_resseqnum1, 520)
        self.assertEqual(graft.orig_insertion1, "")
        self.assertEqual(graft.orig_resseqnum2, None)
        self.assertEqual(graft.orig_insertion2, None)
        self.assertEqual(graft.source_pdbid, "4ijk")
        self.assertEqual(graft.source_chainID, "A")
        self.assertEqual(graft.source_resseqnum1, 1)
        self.assertEqual(graft.source_insertion1, "")
        self.assertEqual(graft.source_resseqnum2, None)
        self.assertEqual(graft.source_insertion2, None)
        self.assertEqual(graft.source_resseqnum3, None)
        self.assertEqual(graft.source_insertion3, None)
        self.assertEqual(graft.obj_id, 0)
        self.assertEqual(graft._yaml_header, 'grafts')
        self.assertEqual(graft._objcat, 'seq')

    def test_graft_from_shortcode(self):
        """
        The shortcode format is "target:source", where:
            - target is in the format "C_RRR[-SSS]"
              - C is the chainID
              - RRR is the residue number and insertion code of the first residue in the target alignment basis
              - SSS is the optional residue number and insertion code of the last residue in the target alignment basis
            - source is in the format "pdbid,C_RRR[#SSS][-TTT]"
              - pdbid is the basename of the PDB file or PDB ID
              - C is the chainID in the source PDB
              - RRR is the residue number and insertion code of the first residue in the source alignment basis
              - SSS is optional, as above
              - TTT is optional, completing RRR-TTT range for entire graft source
        """
        shortcode = "A_520:4ijk,A_1-10"
        graft = Graft.new(shortcode)
        self.assertIsInstance(graft, Graft)
        self.assertEqual(graft.orig_chainID, "A")
        self.assertEqual(graft.orig_resseqnum1, 520)
        self.assertEqual(graft.orig_insertion1, "")
        self.assertEqual(graft.orig_resseqnum2, 520)
        self.assertEqual(graft.orig_insertion2, "")
        self.assertEqual(graft.source_pdbid, "4ijk")
        self.assertEqual(graft.source_chainID, "A")
        self.assertEqual(graft.source_resseqnum1, 1)
        self.assertEqual(graft.source_insertion1, "")
        self.assertEqual(graft.source_resseqnum2, 1)
        self.assertEqual(graft.source_insertion2, "")
        self.assertEqual(graft.source_resseqnum3, 10)
        self.assertEqual(graft.source_insertion3, "")
        self.assertEqual(graft.obj_id, 0)
        # do it again to make sure the counter increments
        graft2 = Graft.new(shortcode)
        self.assertEqual(graft2.obj_id, 1)

class TestGraftList(unittest.TestCase):
    def test_graft_list_creation(self):
        graft_list = GraftList()
        self.assertIsInstance(graft_list, GraftList)
        self.assertEqual(len(graft_list), 0)

    def test_graft_list_addition(self):
        graft = Graft(
            orig_chainID="A",
            orig_resseqnum1=520,
            orig_insertion1="",
            source_pdbid="4ijk",
            source_chainID="A",
            source_resseqnum1=1,
            source_insertion1=""
        )
        graft_list = GraftList()
        graft_list.append(graft)
        self.assertEqual(len(graft_list), 1)
        self.assertEqual(graft_list[0], graft)

    def test_graft_list_assignment(self):
        graft = Graft(
            orig_chainID="A",
            orig_resseqnum1=520,
            orig_insertion1="",
            source_pdbid="4ijk",
            source_chainID="A",
            source_resseqnum1=1,
            source_insertion1=""
        )
        graft_list = GraftList()
        graft_list.append(graft)
        self.assertEqual(graft_list[0].orig_chainID, "A")