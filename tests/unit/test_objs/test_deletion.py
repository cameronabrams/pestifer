import unittest
from pestifer.objs.deletion import Deletion, DeletionList
from pestifer.objs.resid import ResID

class TestDeletion(unittest.TestCase):
    def test_deletion_creation(self):
        deletion = Deletion(
            chainID="A",
            resid1=ResID(10),
            resid2=ResID(20)
        )
        self.assertIsInstance(deletion, Deletion)
        self.assertEqual(deletion.chainID, "A")
        self.assertEqual(deletion.resid1.resid, 10)
        self.assertEqual(deletion.resid2.resid, 20)
        self.assertEqual(deletion._yaml_header, 'deletions')
        self.assertEqual(deletion._objcat, 'seq')
        self.assertEqual(repr(deletion), "Deletion(chainID='A', resid1=ResID(resseqnum=10), resid2=ResID(resseqnum=20))")

    def test_deletion_from_obj(self):
        deletion1 = Deletion(
            chainID="A",
            resid1=ResID(10),
            resid2=ResID(20)
        )
        deletion2 = deletion1.copy()
        self.assertIsInstance(deletion2, Deletion)
        self.assertTrue(deletion1 is not deletion2)
        self.assertEqual(deletion1, deletion2)
        self.assertIsInstance(deletion2, Deletion)
        self.assertEqual(deletion2.chainID, "A")
        self.assertEqual(deletion2.resid1.resid, 10)
        self.assertEqual(deletion2.resid2.resid, 20)

    def test_deletion_from_shortcode(self):
        deletion = Deletion("A:10-20")
        self.assertIsInstance(deletion, Deletion)
        self.assertEqual(deletion.chainID, "A")
        self.assertEqual(deletion.resid1.resid, 10)
        self.assertEqual(deletion.resid2.resid, 20)

    def test_deletion_to_input_string(self):
        deletion = Deletion(
            chainID="A",
            resid1=ResID(10),
            resid2=ResID('20J')
        )
        input_string = deletion.shortcode()
        self.assertEqual(input_string, "A:10-20J")
    
    def test_deletion_list_creation(self):
        deletion_list = DeletionList()
        self.assertIsInstance(deletion_list, DeletionList)
        self.assertEqual(len(deletion_list), 0)