import unittest
from pestifer.objs.substitution import Substitution, SubstitutionList

class TestSubstitution(unittest.TestCase):
    
    def test_substitution_creation(self):
        substitution = Substitution(
            chainID="A",
            resseqnum1=10,
            insertion1="",
            resseqnum2=20,
            insertion2="",
            subseq="G"
        )
        self.assertIsInstance(substitution, Substitution)
        self.assertEqual(substitution.chainID, "A")
        self.assertEqual(substitution.resseqnum1, 10)
        self.assertEqual(substitution.insertion1, "")
        self.assertEqual(substitution.resseqnum2, 20)
        self.assertEqual(substitution.insertion2, "")
        self.assertEqual(substitution.subseq, "G")
        self.assertEqual(substitution._yaml_header, 'substitutions')
        self.assertEqual(substitution._objcat, 'seq')
        self.assertEqual(substitution.describe(), "Substitution(chainID=A, resseqnum1=10, insertion1=, resseqnum2=20, insertion2=, subseq=G)")

    def test_substitution_from_shortcode(self):
        substitution = Substitution.new("A:10-20,G")
        self.assertIsInstance(substitution, Substitution)
        self.assertEqual(substitution.chainID, "A")
        self.assertEqual(substitution.resseqnum1, 10)
        self.assertEqual(substitution.insertion1, "")
        self.assertEqual(substitution.resseqnum2, 20)
        self.assertEqual(substitution.insertion2, "")
        self.assertEqual(substitution.subseq, "G")

    def test_substitution_to_input_string(self):
        substitution = Substitution(
            chainID="A",
            resseqnum1=10,
            insertion1="",
            resseqnum2=20,
            insertion2="J",
            subseq="G"
        )
        input_string = substitution.to_input_string()
        self.assertEqual(input_string, "A:10-20J,G")

    def test_substitution_list_creation(self):
        substitution_list = SubstitutionList()
        self.assertIsInstance(substitution_list, SubstitutionList)
        self.assertEqual(len(substitution_list), 0)