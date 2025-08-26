# Author: Cameron F. Abrams, <cfa22@drexel.edu>
# with substantial contributions from Chat GPT 4o

import unittest
from pestifer.objs.cleavagesite import CleavageSite, CleavageSiteList
from pestifer.objs.resid import ResID

class TestCleavageSite(unittest.TestCase):
    
    def test_cleavage_site_creation(self):
        cleavage_site = CleavageSite(
            chainID="A",
            resid1=ResID(10),
            resid2=ResID(20)
        )
        self.assertIsInstance(cleavage_site, CleavageSite)
        self.assertEqual(cleavage_site.chainID, "A")
        self.assertEqual(cleavage_site.resid1.resid, 10)
        self.assertEqual(cleavage_site.resid2.resid, 20)
        self.assertEqual(cleavage_site._yaml_header, 'cleavages')
        self.assertEqual(cleavage_site._objcat, 'seq')

    def test_cleavage_site_from_shortcode(self):
        cleavage_site = CleavageSite("A:10-20")
        self.assertIsInstance(cleavage_site, CleavageSite)
        self.assertEqual(cleavage_site.chainID, "A")
        self.assertEqual(cleavage_site.resid1.resid, 10)
        self.assertEqual(cleavage_site.resid2.resid, 20)

    def test_cleavage_site_to_shortcode(self):
        cleavage_site = CleavageSite(
            chainID="A",
            resid1=ResID(10),
            resid2=ResID('20J')
        )
        shortcode = cleavage_site.shortcode()
        self.assertEqual(shortcode, "A:10-20J")

        self.assertEqual(repr(cleavage_site), "CleavageSite(chainID='A', resid1=ResID(resseqnum=10), resid2=ResID(resseqnum=20, insertion='J'))")

    def test_cleavage_site_list_creation(self):
        cleavage_list = CleavageSiteList()
        self.assertIsInstance(cleavage_list, CleavageSiteList)
        self.assertEqual(len(cleavage_list), 0)