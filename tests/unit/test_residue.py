import unittest
from pestifer.residue import Residue
from pestifer.config import ConfigSetup
from pestifer.molecule import Molecule

class TestResidue(unittest.TestCase):
    def setUp(self):
        return super().setUp()
    def test_from_atom(self):
        self.config=ConfigSetup('')
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

    def test_segtypes(self):
        self.config=ConfigSetup('')
        m=Molecule.from_rcsb(pdb_code='4tvp')
        p=m.Residues.get(segtype='PROTEIN')
        pc=[]
        for r in p:
            if not r.chainID in pc:
                pc.append(r.chainID)
        self.assertEqual(pc,['G', 'B', 'L', 'H', 'D', 'E'])
        i=m.Residues.get(segtype='ION')
        ic=[]
        for r in i:
            if not r.chainID in ic:
                ic.append(r.chainID)
        self.assertEqual(ic,['G','B','L'])
        w=m.Residues.get(segtype='WATER')
        wc=[]
        for r in w:
            if not r.chainID in wc:
                wc.append(r.chainID)
        self.assertEqual(wc,['G', 'B', 'L', 'H', 'D', 'E'])
        g=m.Residues.get(segtype='GLYCAN')
        gc=[]
        for r in g:
            if not r.chainID in gc:
                gc.append(r.chainID)
        self.assertEqual(gc,['A','C','F','I','J','K','M','N','O','P','Q','R','S','T','G','B','H'])
        l=m.Residues.get(segtype='LIGAND')
        lc=[]
        for r in l:
            if not r.chainID in lc:
                lc.append(r.chainID)
        self.assertEqual(lc,[])
        o=m.Residues.get(segtype='OTHER')
        oc=[]
        for r in o:
            if not r.chainID in oc:
                oc.append(r.chainID)
        self.assertEqual(oc,[])
    def test_residuelist(self):
        self.config=ConfigSetup('')
        m=Molecule.from_rcsb(pdb_code='6m0j')
        r=m.Residues.get_residue(resseqnum=427,chainID='A')
        self.assertEqual(r.name,'ASP')
        a=m.Residues.get_atom(atname='OD2',resseqnum=427,chainID='A')
        self.assertEqual(a.serial,3302)
