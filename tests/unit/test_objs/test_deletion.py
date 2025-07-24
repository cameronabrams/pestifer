import unittest
from pestifer.objs.deletion import Deletion, DeletionList

class TestDeletion(unittest.TestCase):
    def test_deletion_creation(self):
        deletion = Deletion(
            chainID="A",
            resseqnum1=10,
            insertion1="",
            resseqnum2=20,
            insertion2=""
        )
        self.assertIsInstance(deletion, Deletion)
        self.assertEqual(deletion.chainID, "A")
        self.assertEqual(deletion.resseqnum1, 10)
        self.assertEqual(deletion.insertion1, "")
        self.assertEqual(deletion.resseqnum2, 20)
        self.assertEqual(deletion.insertion2, "")
        self.assertEqual(deletion._yaml_header, 'deletions')
        self.assertEqual(deletion._objcat, 'seq')
        self.assertEqual(deletion.describe(), "Deleteion(chainID=A, resseqnum1=10, insertion1=, resseqnum2=20, insertion2=)")

    def test_deletion_from_obj(self):
        deletion1 = Deletion(
            chainID="A",
            resseqnum1=10,
            insertion1="",
            resseqnum2=20,
            insertion2=""
        )
        deletion2 = Deletion.from_input(deletion1)
        self.assertTrue(deletion1 is not deletion2)
        self.assertEqual(deletion1, deletion2)
        self.assertIsInstance(deletion2, Deletion)
        self.assertEqual(deletion2.chainID, "A")
        self.assertEqual(deletion2.resseqnum1, 10)
        self.assertEqual(deletion2.insertion1, "")
        self.assertEqual(deletion2.resseqnum2, 20)
        self.assertEqual(deletion2.insertion2, "")

    def test_deletion_from_shortcode(self):
        shortcode = Deletion.Adapter.from_string("A:10-20")
        deletion = Deletion.from_input(shortcode)
        self.assertIsInstance(deletion, Deletion)
        self.assertEqual(deletion.chainID, "A")
        self.assertEqual(deletion.resseqnum1, 10)
        self.assertEqual(deletion.insertion1, "")
        self.assertEqual(deletion.resseqnum2, 20)
        self.assertEqual(deletion.insertion2, "")

    def test_deletion_to_input_string(self):
        deletion = Deletion(
            chainID="A",
            resseqnum1=10,
            insertion1="",
            resseqnum2=20,
            insertion2="J"
        )
        input_string = deletion.to_input_string()
        self.assertEqual(input_string, "A:10-20J")
    
    def test_deletion_list_creation(self):
        deletion_list = DeletionList()
        self.assertIsInstance(deletion_list, DeletionList)
        self.assertEqual(len(deletion_list), 0)