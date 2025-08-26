import unittest
from pestifer.objs.insertion import Insertion, InsertionList
from pestifer.objs.resid import ResID

class TestInsertion(unittest.TestCase):
    def test_insertion_creation(self):
        insertion = Insertion(
            chainID="A",
            resid=ResID(10),
            sequence="GGG",
            integer_increment=False
        )
        self.assertIsInstance(insertion, Insertion)
        self.assertEqual(insertion.chainID, "A")
        self.assertEqual(insertion.resid.resid, 10)
        self.assertEqual(insertion.sequence, "GGG")
        self.assertFalse(insertion.integer_increment)
        self.assertEqual(insertion._yaml_header, 'insertions')
        self.assertEqual(insertion._objcat, 'seq')
        self.assertEqual(repr(insertion), "Insertion(chainID='A', resid=ResID(resseqnum=10), sequence='GGG', integer_increment=False)")

    def test_insertion_from_shortcode(self):
        insertion = Insertion("A,10,GGG")
        self.assertIsInstance(insertion, Insertion)
        self.assertEqual(insertion.chainID, "A")
        self.assertEqual(insertion.resid.resid, 10)
        self.assertEqual(insertion.sequence, "GGG")
        self.assertFalse(insertion.integer_increment)

    def test_insertion_from_shortcode_with_increment(self):
        insertion = Insertion("A,10+,GRETA")
        self.assertIsInstance(insertion, Insertion)
        self.assertEqual(insertion.chainID, "A")
        self.assertEqual(insertion.resid.resid, 10)
        self.assertEqual(insertion.sequence, "GRETA")
        self.assertTrue(insertion.integer_increment)


class TestInsertionList(unittest.TestCase):
    def test_insertion_list_creation(self):
        insertion_list = InsertionList()
        self.assertIsInstance(insertion_list, InsertionList)
        self.assertEqual(len(insertion_list), 0)