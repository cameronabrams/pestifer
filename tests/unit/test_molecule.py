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



    def test_links(self):
        self.config=ConfigSetup('')
        m=Molecule.from_rcsb(pdb_code='4zmj')
        l=m.Links[0]
        self.assertEqual(l.residue1.segtype,'PROTEIN')
        self.assertEqual(l.residue1.name,'ASN')
        self.assertEqual(l.residue1.chainID,'G')
        self.assertEqual(l.residue1.resseqnum,156)
        self.assertEqual(l.atom1.name,'ND2')
        self.assertEqual(l.atom1.altloc,'')
        self.assertEqual(l.residue2.segtype,'GLYCAN')
        self.assertEqual(l.residue2.name,'NAG')
        self.assertEqual(l.residue2.chainID,'G')
        self.assertEqual(l.residue2.resseqnum,615)
        self.assertEqual(l.atom2.name,'C1')
        self.assertEqual(l.atom2.altloc,'')
        self.assertTrue(l.residue2 in l.residue1.down)
        self.assertTrue(l.residue1 in l.residue2.up)
        self.assertTrue(l in l.residue1.downlink)
        self.assertTrue(l in l.residue2.uplink)
        l=m.Links[-1]
        self.assertEqual(l.residue1.segtype,'GLYCAN')
        self.assertEqual(l.residue1.name,'NAG')
        self.assertEqual(l.residue1.chainID,'D')
        self.assertEqual(l.residue1.resseqnum,1)
        self.assertEqual(l.atom1.name,'O4')
        self.assertEqual(l.atom1.altloc,'')
        self.assertEqual(l.residue2.segtype,'GLYCAN')
        self.assertEqual(l.residue2.name,'NAG')
        self.assertEqual(l.residue2.chainID,'D')
        self.assertEqual(l.residue2.resseqnum,2)
        self.assertEqual(l.atom2.name,'C1')
        self.assertEqual(l.atom2.altloc,'')

    def test_bioassemb(self):
        self.config=ConfigSetup('')
        m=Molecule.from_rcsb(pdb_code='4zmj')
        self.assertEqual(2,len(m.BioAssemb))