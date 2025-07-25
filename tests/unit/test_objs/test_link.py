import unittest
from pestifer.objs.link import Link, LinkList

class TestLink(unittest.TestCase):

    def test_link_creation(self):
        link = Link(
            chainID1="A",
            resseqnum1=1,
            insertion1="",
            name1="N",
            chainID2="B",
            resseqnum2=2,
            insertion2="",
            name2="CA"
        )
        self.assertIsInstance(link, Link)
        self.assertEqual(link.chainID1, "A")
        self.assertEqual(link.resseqnum1, 1)
        self.assertEqual(link.insertion1, "")
        self.assertEqual(link.name1, "N")
        self.assertEqual(link.chainID2, "B")
        self.assertEqual(link.resseqnum2, 2)
        self.assertEqual(link.insertion2, "")
        self.assertEqual(link.name2, "CA")
        self.assertEqual(link._yaml_header, 'links')
        self.assertEqual(link._objcat, 'topol')
        self.assertEqual(link.describe(), "Link between N and CA (A:1-B:2)")

    def test_link_from_shortcode(self):
        link = Link.new("A_1_N-B_2_CA")
        self.assertIsInstance(link, Link)
        self.assertEqual(link.chainID1, "A")
        self.assertEqual(link.resseqnum1, 1)
        self.assertEqual(link.insertion1, "")
        self.assertEqual(link.name1, "N")
        self.assertEqual(link.chainID2, "B")
        self.assertEqual(link.resseqnum2, 2)
        self.assertEqual(link.insertion2, "")
        self.assertEqual(link.name2, "CA")

