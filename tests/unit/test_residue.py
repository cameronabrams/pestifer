import unittest
from pestifer.residue import Residue
from pestifer.mods import Atom, Missing

class TestResidue(unittest.TestCase):
    def setUp(self):
        self.pdbrecord='REMARK 465     GLU G   185A  '
        pr1='ATOM    981  N   GLU G 164A     10.462 130.792 -22.224  1.00 70.92           N  '
        pr2='ATOM    984  O   GLU G 164A     10.462 130.792 -22.224  1.00 70.92           N  '
        pr3='ATOM    984  O   GLU G 164      10.462 130.792 -22.224  1.00 70.92           N  '
        pr4='ATOM    984  O   GLU G 165      10.462 130.792 -22.224  1.00 70.92           N  '
        pr5='ATOM    984  O   GLU G 165A     10.462 130.792 -22.224  1.00 70.92           N  '
        self.missing=Missing.from_pdbrecord(self.pdbrecord)
        self.atom1=Atom.from_pdbrecord(pdbrecord=pr1)
        self.atom2=Atom.from_pdbrecord(pdbrecord=pr2)
        self.atom3=Atom.from_pdbrecord(pdbrecord=pr3)
        self.atom4=Atom.from_pdbrecord(pdbrecord=pr4)
        self.atom5=Atom.from_pdbrecord(pdbrecord=pr5)
        return super().setUp()
    def test_from_atom(self):
        res=Residue.from_atom(self.atom1)
        self.assertEqual(res.resseqnum,self.atom1.resseqnum)
        self.assertEqual(len(res.atoms),1)
        res.add_atom(self.atom2)
        self.assertEqual(len(res.atoms),2)
        res.add_atom(self.atom3)
        self.assertEqual(len(res.atoms),2)
    def test_ordinality(self):
        res1=Residue.from_atom(self.atom3)
        res2=Residue.from_atom(self.atom4)
        res3=Residue.from_atom(self.atom5)
        self.assertTrue(res1<res2)
        self.assertTrue(res2<res3)
    def test_set_chainid(self):
        res=Residue.from_atom(self.atom1)
        res.add_atom(self.atom2)
        res.set_chainID('X')
        self.assertEqual(res.chainID,'X')
        for i in range(len(res.atoms)):
            self.assertEqual(res.atoms[i].chainID,'X')