import unittest
from pestifer.residue import Residue
from pestifer.config import Config
from pestifer.molecule import Molecule

class TestResidue(unittest.TestCase):
    def test_from_atom(self):
        m=Molecule(source='1gc1',config=Config())
        au=m.asymmetric_unit
        residues=[]
        a=au.Atoms[0]
        r=Residue(a)
        self.assertEqual(a.chainID,r.chainID)
        residues.append(r)
        for a in au.Atoms[1:]:
            if not any([x.add_atom(a) for x in residues]):
                r=Residue(a)
                residues.append(r)
        self.assertEqual(len(residues),1538)

    def test_segtypes(self):
        m=Molecule(source='4tvp',config=Config())
        au=m.asymmetric_unit
        p=au.Residues.get(segtype='PROTEIN')
        pc=[]
        for r in p:
            if not r.chainID in pc:
                pc.append(r.chainID)
        self.assertEqual(pc,['G', 'B', 'L', 'H', 'D', 'E'])
        i=au.Residues.get(segtype='ION')
        ic=[]
        for r in i:
            if not r.chainID in ic:
                ic.append(r.chainID)
        self.assertEqual(ic,['G','B','L'])
        w=au.Residues.get(segtype='WATER')
        wc=[]
        for r in w:
            if not r.chainID in wc:
                wc.append(r.chainID)
        self.assertEqual(wc,['G', 'B', 'L', 'H', 'D', 'E'])
        g=au.Residues.get(segtype='GLYCAN')
        gc=[]
        for r in g:
            if not r.chainID in gc:
                gc.append(r.chainID)
        self.assertEqual(gc,['A','C','F','I','J','K','M','N','O','P','Q','R','S','T','G','B','H'])
        l=au.Residues.get(segtype='LIGAND')
        lc=[]
        for r in l:
            if not r.chainID in lc:
                lc.append(r.chainID)
        self.assertEqual(lc,[])
        o=au.Residues.get(segtype='OTHER')
        oc=[]
        for r in o:
            if not r.chainID in oc:
                oc.append(r.chainID)
        self.assertEqual(oc,[])

    def test_residuelist(self):
        m=Molecule(source='6m0j',config=Config())
        au=m.asymmetric_unit
        r=au.Residues.get_residue(resseqnum=427,chainID='A')
        self.assertEqual(r.name,'ASP')
        a=au.Residues.get_atom('OD2',resseqnum=427,chainID='A')
        self.assertEqual(a.serial,3302)
