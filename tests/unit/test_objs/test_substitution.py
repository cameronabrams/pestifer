import unittest
from pestifer.objs.substitution import Substitution, SubstitutionList
from pestifer.objs.resid import ResID, ResIDList

class TestSubstitution(unittest.TestCase):
    
    def test_substitution_creation(self):
        substitution = Substitution(
            chainID="A",
            resid1=ResID(resseqnum=10),
            resid2=ResID(resseqnum=20),
            subseq="G"
        )
        self.assertIsInstance(substitution, Substitution)
        self.assertEqual(substitution.chainID, "A")
        self.assertEqual(substitution.resid1.resseqnum, 10)
        self.assertEqual(substitution.resid1.insertion, None)
        self.assertEqual(substitution.resid2.resseqnum, 20)
        self.assertEqual(substitution.resid2.insertion, None)
        self.assertEqual(substitution.subseq, "G")
        self.assertEqual(substitution._yaml_header, 'substitutions')
        self.assertEqual(substitution._objcat, 'seq')
        self.assertEqual(repr(substitution), "Substitution(chainID='A', resid1=ResID(resseqnum=10), resid2=ResID(resseqnum=20), subseq='G')")

    def test_substitution_from_shortcode(self):
        substitution = Substitution("A:10-20,G")
        self.assertIsInstance(substitution, Substitution)
        self.assertEqual(substitution.chainID, "A")
        self.assertEqual(substitution.resid1, ResID(10))
        self.assertEqual(substitution.resid2, ResID(20))
        self.assertEqual(substitution.subseq, "G")

    def test_substitution_to_input_string(self):
        substitution = Substitution(
            chainID="A",
            resid1=ResID(resseqnum=10),
            resid2=ResID(resseqnum=20, insertion="J"),
            subseq="G"
        )
        input_string = substitution.shortcode()
        self.assertEqual(input_string, "A:10-20J,G")

    def test_substitution_list_creation(self):
        substitution_list = SubstitutionList()
        self.assertIsInstance(substitution_list, SubstitutionList)
        self.assertEqual(len(substitution_list), 0)