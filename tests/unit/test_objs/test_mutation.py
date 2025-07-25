import unittest
from pestifer.objs.mutation import Mutation, MutationList

class TestMutation(unittest.TestCase):
    def test_mutation_creation(self):
        mutation = Mutation(
            chainID="A",
            origresname="ALA",
            resseqnum=10,
            insertion="",
            newresname="GLY",
            typekey="user",
            pdbx_auth_seq_num=10
        )
        self.assertIsInstance(mutation, Mutation)
        self.assertEqual(mutation.chainID, "A")
        self.assertEqual(mutation.origresname, "ALA")
        self.assertEqual(mutation.resseqnum, 10)
        self.assertEqual(mutation.insertion, "")
        self.assertEqual(mutation.newresname, "GLY")
        self.assertEqual(mutation.typekey, "user")
        self.assertEqual(mutation.pdbx_auth_seq_num, 10)
        self.assertEqual(mutation._yaml_header, 'mutations')
        self.assertEqual(mutation._objcat, 'seq')
        self.assertEqual(mutation.describe(), "Mutation(chainID=A, origresname=ALA, resseqnum=10, insertion=, newresname=GLY, typekey=user, pdbx_auth_seq_num=10)")

    def test_mutation_from_shortcode(self):
        mutation = Mutation.new("C:ALA,10,GLY")
        self.assertIsInstance(mutation, Mutation)
        self.assertEqual(mutation.chainID, "C")
        self.assertEqual(mutation.origresname, "ALA")
        self.assertEqual(mutation.resseqnum, 10)
        self.assertEqual(mutation.insertion, "")
        self.assertEqual(mutation.newresname, "GLY")
        self.assertEqual(mutation.typekey, "user")
        self.assertEqual(mutation.pdbx_auth_seq_num, None)