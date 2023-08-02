import unittest
import os
from pestifer.chain import Chain
from pestifer.molecule import Molecule
from pestifer.residue import Residue
from pestifer.config import Config
class TestChain(unittest.TestCase):

    def test1(self):
        c=Chain({'chainID':'A'})
        self.assertEqual(c.chainID,'A')
        self.assertTrue(hasattr(c,'residues'))

    def test_from_residue(self):
        m=Molecule.from_rcsb(pdb_code='1gc1')
        r=m.Residues[0]
        c=Chain.from_residue(r,parent_molecule=m)
        self.assertEqual(r.chainID,c.chainID)
        chains=[]
        chains.append(c)
        for r in m.Residues[1:]:
            if not any([x.add_residue(r) for x in chains]):
                c=Chain.from_residue(r,parent_molecule=m)
                chains.append(c)
        self.assertEqual(len(chains),8)
        expected=['G', 'C', 'L', 'H', 'A', 'B', 'D', 'E']
        chainIDs=[x.chainID for x in chains]
        self.assertEqual(expected,chainIDs)
        c=chains[0]
        self.assertEqual(c.get_molid(),m.molid)
        c.sort_residues()
        n=min([x.resseqnum for x in c.residues])
        self.assertEqual(n,c.residues[0].resseqnum)
        # residues.append(r)
        # for a in m.Atoms[1:]:
        #     if not r.add_atom(a):
        #         r=Residue.from_atom(a)
        #         residues.append(r)
        # self.assertEqual(len(residues),1538)