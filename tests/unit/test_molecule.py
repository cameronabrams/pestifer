import unittest
import os
from pestifer.molecule import Molecule
from pestifer.command import Command
class TestMolecule(unittest.TestCase):
    def setUp(self):
        return super().setUp()
    def test_molecule(self):
        m=Molecule.from_pdb(pdb_code='1gc1')
        self.assertEqual(m.pdb_code,'1gc1')
        c=Command('wc -l 1gc1.pdb')
        c.run()
        nl=int(c.stdout.split()[0])
        self.assertEqual(nl,len(m.pdb_parser.pdb_lines))
        os.remove('1gc1.pdb')