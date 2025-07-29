import unittest
from pestifer.objs.link import Link, LinkList
from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict, PDBRecordList
from pestifer.util.cifutil import CIFload, CIFdict
from pathlib import Path
from mmcif.api.PdbxContainers import DataContainer
from mmcif.api.DataCategory import DataCategory

import logging
logger = logging.getLogger(__name__)
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
        self.assertEqual(link.patchhead,1)
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

class TestLinkList(unittest.TestCase):

    def test_link_list_from_pdb(self):
        p = PDBParser(filepath='data/4zmj.pdb').parse().parsed
        self.assertIsInstance(p, PDBRecordDict)
        L = LinkList.from_pdb(p)
        self.assertIsInstance(L, LinkList)
        self.assertGreater(len(L), 0)

    def test_link_list_from_cif(self):
        cif_data = CIFload(Path('data/4zmj.cif'))
        self.assertIsInstance(cif_data, DataContainer)
        L = LinkList.from_cif(cif_data)
        self.assertIsInstance(L, LinkList)
        self.assertGreater(len(L), 0)

