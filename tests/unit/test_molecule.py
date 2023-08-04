import unittest
from pestifer.molecule import Molecule
from pestifer.config import ConfigSetup
class TestMolecule(unittest.TestCase):
    def setUp(self):
        return super().setUp()
    def test_molecule(self):
        self.config=ConfigSetup('')
        m=Molecule.from_rcsb(pdb_code='1gc1')
        self.assertEqual(m.pdb_code,'1gc1')
        self.assertEqual(len(m.Atoms),7877)
        self.assertEqual(len(m.SSBonds),14)
        self.assertEqual(len(m.Mutations),2)
        self.assertEqual(len(m.Conflicts),1)
        self.assertEqual(len(m.Missing),28)
        self.assertEqual(len(m.Links),15)
        self.assertEqual(len(m.Residues),1566)
        self.assertEqual(m.Residues[0].segtype,'PROTEIN')
        waterres=[r for r in m.Residues if r.segtype=='WATER']
        self.assertEqual(len(waterres),603)
        glycanres=[r for r in m.Residues if r.segtype=='GLYCAN']
        self.assertEqual(len(glycanres),15)
        g0=glycanres[0]
        self.assertEqual(g0.name,'NAG')

    def test_bioassemb(self):
        self.config=ConfigSetup('')
        m=Molecule.from_rcsb(pdb_code='4zmj')
        self.assertEqual(2,len(m.BioAssemb))