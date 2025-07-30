# Author: Cameron F. Abrams, <cfa22@drexel.edu>
# with substantial contributions from Chat GPT 4o

import unittest
from pestifer.objs.cleavagesite import CleavageSite, CleavageSiteList

class TestCleavageSite(unittest.TestCase):
    
    def test_cleavage_site_creation(self):
        cleavage_site = CleavageSite(
            chainID="A",
            resseqnum1=10,
            insertion1="",
            resseqnum2=20,
            insertion2=""
        )
        self.assertIsInstance(cleavage_site, CleavageSite)
        self.assertEqual(cleavage_site.chainID, "A")
        self.assertEqual(cleavage_site.resseqnum1, 10)
        self.assertEqual(cleavage_site.insertion1, "")
        self.assertEqual(cleavage_site.resseqnum2, 20)
        self.assertEqual(cleavage_site.insertion2, "")
        self.assertEqual(cleavage_site._yaml_header, 'cleavages')
        self.assertEqual(cleavage_site._objcat, 'seq')

    def test_cleavage_site_from_shortcode(self):
        cleavage_site = CleavageSite("A:10-20")
        self.assertIsInstance(cleavage_site, CleavageSite)
        self.assertEqual(cleavage_site.chainID, "A")
        self.assertEqual(cleavage_site.resseqnum1, 10)
        self.assertEqual(cleavage_site.insertion1, "")
        self.assertEqual(cleavage_site.resseqnum2, 20)
        self.assertEqual(cleavage_site.insertion2, "")

    def test_cleavage_site_to_shortcode(self):
        cleavage_site = CleavageSite(
            chainID="A",
            resseqnum1=10,
            insertion1="",
            resseqnum2=20,
            insertion2="J"
        )
        shortcode = cleavage_site.shortcode()
        self.assertEqual(shortcode, "A:10-20J")

        self.assertEqual(repr(cleavage_site), "CleavageSite(chainID='A', resseqnum1=10, insertion1='', resseqnum2=20, insertion2='J')")

    def test_cleavage_site_list_creation(self):
        cleavage_list = CleavageSiteList()
        self.assertIsInstance(cleavage_list, CleavageSiteList)
        self.assertEqual(len(cleavage_list), 0)