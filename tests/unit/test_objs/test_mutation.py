import unittest
from pestifer.objs.mutation import Mutation, MutationList
from pidibble.pdbparse import PDBParser
from pestifer.molecule.residue import Residue
from pestifer.objs.seqadv import Seqadv
from pestifer.objs.resid import ResID

from pathlib import Path

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

        mutation2 = Mutation("C:ALA10GLY")
        self.assertEqual(mutation2.chainID, "C")
        self.assertEqual(mutation2.origresname, "ALA")
        self.assertEqual(mutation2.resid, ResID(10))
        self.assertEqual(mutation2.newresname, "GLY")
        self.assertEqual(mutation2.typekey, "user")
        self.assertEqual(mutation2.pdbx_auth_seq_num, None)

        mutation3 = Mutation("C:A,10,G")
        self.assertEqual(mutation3.chainID, "C")
        self.assertEqual(mutation3.origresname, "ALA")
        self.assertEqual(mutation3.resid, ResID(10))
        self.assertEqual(mutation3.newresname, "GLY")
        self.assertEqual(mutation3.typekey, "user")
        self.assertEqual(mutation3.pdbx_auth_seq_num, None)

        mutation4 = Mutation("C:A10G")
        self.assertEqual(mutation4.chainID, "C")
        self.assertEqual(mutation4.origresname, "ALA")
        self.assertEqual(mutation4.resid, ResID(10))
        self.assertEqual(mutation4.newresname, "GLY")
        self.assertEqual(mutation4.typekey, "user")
        self.assertEqual(mutation4.pdbx_auth_seq_num, None)

    def test_mutation_from_seqadv(self):
        inputs_dir = Path(__file__).parents[2] / "inputs"
        p=PDBParser(filepath=inputs_dir / '4zmj.pdb').parse()
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