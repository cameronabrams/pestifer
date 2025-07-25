import unittest
from pestifer.objs.insertion import Insertion, InsertionList

class TestInsertion(unittest.TestCase):
    def test_insertion_creation(self):
        insertion = Insertion(
            chainID="A",
            resseqnum=10,
            insertion="",
            sequence="GGG",
            integer_increment=False
        )
        self.assertIsInstance(insertion, Insertion)
        self.assertEqual(insertion.chainID, "A")
        self.assertEqual(insertion.resseqnum, 10)
        self.assertEqual(insertion.insertion, "")
        self.assertEqual(insertion.sequence, "GGG")
        self.assertFalse(insertion.integer_increment)
        self.assertEqual(insertion._yaml_header, 'insertions')
        self.assertEqual(insertion._objcat, 'seq')
        self.assertEqual(insertion.describe(), "Insertion(chainID=A, resseqnum=10, insertion=, sequence=GGG, integer_increment=False)")

    def test_insertion_from_shortcode(self):
        insertion = Insertion.new("A,10,GGG")
        self.assertIsInstance(insertion, Insertion)
        self.assertEqual(insertion.chainID, "A")
        self.assertEqual(insertion.resseqnum, 10)
        self.assertEqual(insertion.insertion, "")
        self.assertEqual(insertion.sequence, "GGG")
        self.assertFalse(insertion.integer_increment)

    def test_insertion_from_shortcode_with_increment(self):
        insertion = Insertion.new("A,10+,GRETA")
        self.assertIsInstance(insertion, Insertion)
        self.assertEqual(insertion.chainID, "A")
        self.assertEqual(insertion.resseqnum, 10)
        self.assertEqual(insertion.insertion, "")
        self.assertEqual(insertion.sequence, "GRETA")
        self.assertTrue(insertion.integer_increment)