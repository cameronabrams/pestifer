import unittest
from pestifer.objs.mutation import Mutation, MutationList
from pidibble.pdbparse import PDBParser
from pestifer.objs.seqadv import Seqadv
from pestifer.objs.resid import ResID

class TestMutation(unittest.TestCase):
    def test_mutation_creation(self):
        mutation = Mutation(
            chainID="A",
            origresname="ALA",
            resid=ResID(10),
            newresname="GLY",
            typekey="user",
            pdbx_auth_seq_num=10
        )
        self.assertIsInstance(mutation, Mutation)
        self.assertEqual(mutation.chainID, "A")
        self.assertEqual(mutation.origresname, "ALA")
        self.assertEqual(mutation.resid, ResID(10))
        self.assertEqual(mutation.newresname, "GLY")
        self.assertEqual(mutation.typekey, "user")
        self.assertEqual(mutation.pdbx_auth_seq_num, 10)
        self.assertEqual(mutation._yaml_header, 'mutations')
        self.assertEqual(mutation._objcat, 'seq')
        self.assertEqual(repr(mutation), "Mutation(chainID='A', origresname='ALA', resid=ResID(resseqnum=10), newresname='GLY', typekey='user', pdbx_auth_seq_num=10)")

    def test_mutation_from_shortcode(self):
        mutation = Mutation("C:ALA,10,GLY")
        self.assertIsInstance(mutation, Mutation)
        self.assertEqual(mutation.chainID, "C")
        self.assertEqual(mutation.origresname, "ALA")
        self.assertEqual(mutation.resid, ResID(10))
        self.assertEqual(mutation.newresname, "GLY")
        self.assertEqual(mutation.typekey, "user")
        self.assertEqual(mutation.pdbx_auth_seq_num, None)

    def test_mutation_from_seqadv(self):
        p=PDBParser(filepath='data/4zmj.pdb').parse()
        pr=p.parsed[Seqadv._PDB_keyword][0]
        seqadv = Seqadv(pr)
        self.assertIsInstance(seqadv, Seqadv)
        m=Mutation(seqadv)
        self.assertIsInstance(m, Mutation)
        self.assertEqual(m.chainID, seqadv.chainID)
        self.assertEqual(m.origresname, seqadv.resname)
        self.assertEqual(m.resid, seqadv.resid)
        self.assertEqual(m.newresname, seqadv.dbRes)
        self.assertEqual(m.typekey, seqadv.typekey)
        self.assertEqual(m.pdbx_auth_seq_num, seqadv.pdbx_auth_seq_num)