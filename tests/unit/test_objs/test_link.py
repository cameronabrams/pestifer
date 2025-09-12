# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import unittest

from mmcif.api.PdbxContainers import DataContainer

from pathlib import Path

from pidibble.pdbparse import PDBParser
from pidibble.pdbrecord import PDBRecordDict

from pestifer.molecule.atom import AtomList
from pestifer.molecule.residue import Residue, ResidueList
from pestifer.objs.link import Link, LinkList
from pestifer.objs.resid import ResID
from pestifer.psfutil.psfpatch import PSFLinkPatch
from pestifer.util.cifutil import CIFload

logger = logging.getLogger(__name__)

class TestLink(unittest.TestCase):

    def test_link_creation(self):
        link = Link(
            chainID1="A",
            resid1=ResID(1),
            name1="N",
            chainID2="B",
            resid2=ResID(2),
            name2="CA"
        )
        self.assertIsInstance(link, Link)
        self.assertEqual(link.chainID1, "A")
        self.assertEqual(link.resid1, ResID(1))
        self.assertEqual(link.name1, "N")
        self.assertEqual(link.chainID2, "B")
        self.assertEqual(link.resid2, ResID(2))
        self.assertEqual(link.name2, "CA")
        self.assertEqual(link._yaml_header, 'links')
        self.assertEqual(link._objcat, 'topol')
        # self.assertEqual(link.patchhead,1)
        self.assertEqual(repr(link), "Link(chainID1='A', resid1=ResID(resseqnum=1), name1='N', chainID2='B', resid2=ResID(resseqnum=2), name2='CA')")

    def test_link_from_shortcode(self):
        link = Link("A_1_N-B_2_CA")
        self.assertIsInstance(link, Link)
        self.assertEqual(link.chainID1, "A")
        self.assertEqual(link.resid1, ResID(1))
        self.assertEqual(link.name1, "N")
        self.assertEqual(link.chainID2, "B")
        self.assertEqual(link.resid2, ResID(2))
        self.assertEqual(link.name2, "CA")

    def test_link_from_psflinkpatch(self):
        patch = PSFLinkPatch(['NGLB', 'A:1', 'B:2'])
        link = Link(patch)
        self.assertIsInstance(link, Link)
        self.assertEqual(link.chainID1, 'A')
        self.assertEqual(link.resid1, ResID(1))
        self.assertEqual(link.name1, 'ND2')
        self.assertEqual(link.chainID2, 'B')
        self.assertEqual(link.resid2, ResID(2))
        self.assertEqual(link.name2, 'C1')

class TestLinkList(unittest.TestCase):

    def setUp(self):
        # set the inputs path to PACKAGE/tests/inputs
        self.inputs_path = Path(__file__).parent.parent.parent / 'inputs'

    def test_link_list_from_pdb(self):
        p = PDBParser(filepath=self.inputs_path / "4zmj.pdb").parse().parsed
        self.assertIsInstance(p, PDBRecordDict)
        L = LinkList.from_pdb(p)
        self.assertIsInstance(L, LinkList)
        self.assertGreater(len(L), 0)

    def test_link_list_from_cif(self):
        cif_data = CIFload(self.inputs_path / '4zmj.cif')
        self.assertIsInstance(cif_data, DataContainer)
        L = LinkList.from_cif(cif_data)
        self.assertIsInstance(L, LinkList)
        self.assertGreater(len(L), 0)
            
    def test_link_list_assign_residues(self):
        # make a list of 10 mock residues
        residues = ResidueList([Residue(chainID='A' if i%2==1 else 'B', segname='A' if i%2==1 else 'B', resid=ResID(i), resname='ALA', atoms=AtomList([]),segtype='protein', resolved=True) for i in range(1, 11)])

        # create a LinkList with 6 links
        links = LinkList([
            Link(chainID1='A', resid1=ResID(1), name1='N', chainID2='B', resid2=ResID(2), name2='CA'),
            Link(chainID1='A', resid1=ResID(3), name1='C', chainID2='B', resid2=ResID(4), name2='O'),
            Link(chainID1='A', resid1=ResID(5), name1='CA', chainID2='B', resid2=ResID(6), name2='CB'),
            Link(chainID1='A', resid1=ResID(7), name1='N', chainID2='B', resid2=ResID(8), name2='CA'),
            Link(chainID1='A', resid1=ResID(9), name1='C', chainID2='B', resid2=ResID(10), name2='O'),
            Link(chainID1='A', resid1=ResID(11), name1='CA', chainID2='B', resid2=ResID(12), name2='CB')
        ])

        # assign residues to the links
        ignored_by_ptnr1 = links.assign_objs_to_attr('residue1', residues, chainID='chainID1', resid='resid1')
        ignored_by_ptnr2 = links.assign_objs_to_attr('residue2', residues, chainID='chainID2', resid='resid2')

        self.assertEqual(links[0].residue1, residues[0])
        self.assertEqual(links[0].residue2, residues[1])
        self.assertEqual(links[1].residue1, residues[2])
        self.assertEqual(links[1].residue2, residues[3])
        self.assertEqual(links[2].residue1, residues[4])
        self.assertEqual(links[2].residue2, residues[5])
        self.assertEqual(len(ignored_by_ptnr1), 1)  # last link references residues that are not in the list
        self.assertEqual(len(ignored_by_ptnr2), 0)  # all links have valid residue2 since the previous link was removed by assign_objs_to_attr