import unittest
from pestifer.residue import Residue
from pestifer.mods import Atom, Missing
from pestifer.molecule import Molecule

class TestResidue(unittest.TestCase):
    def setUp(self):
        return super().setUp()
    def test_from_atom(self):
        m=Molecule.from_rcsb(pdb_code='1gc1')
        residues=[]
        a=m.Atoms[0]
        r=Residue.from_atom(a)
        self.assertEqual(a.chainID,r.chainID)
        residues.append(r)
        for a in m.Atoms[1:]:
            if not any([x.add_atom(a) for x in residues]):
                r=Residue.from_atom(a)
                residues.append(r)
        self.assertEqual(len(residues),1538)

    def test_residuelist(self):
        m=Molecule.from_rcsb(pdb_code='6m0j')
        r=m.Residues.get_residue(resseqnum=427,chainID='A')
        self.assertEqual(r.name,'ASP')
        a=m.Residues.get_atom(atname='OD2',resseqnum=427,chainID='A')
        self.assertEqual(a.serial,3302)
