import unittest
from pestifer.objs.mutation import Mutation, MutationList
from pidibble.pdbparse import PDBParser
from pestifer.objs.seqadv import Seqadv

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

    def test_mutation_from_seqadv(self):
        p=PDBParser(filepath='data/4zmj.pdb').parse()
        pr=p.parsed[Seqadv._PDB_keyword][0]
        seqadv = Seqadv.new(pr)
        self.assertIsInstance(seqadv, Seqadv)
        m=Mutation.new(seqadv)
        self.assertIsInstance(m, Mutation)
        self.assertEqual(m.chainID, seqadv.chainID)
        self.assertEqual(m.origresname, seqadv.resname)
        self.assertEqual(m.resseqnum, seqadv.resseqnum)
        self.assertEqual(m.insertion, seqadv.insertion)
        self.assertEqual(m.newresname, seqadv.dbRes)
        self.assertEqual(m.typekey, seqadv.typekey)
        self.assertEqual(m.pdbx_auth_seq_num, seqadv.pdbx_auth_seq_num)