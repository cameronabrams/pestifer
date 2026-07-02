import unittest
from pestifer.psfutil.psfring import ring_check
import os
import logging
import pytest
logger=logging.getLogger(__name__)

class TestRingCheck(unittest.TestCase):
    def test_ring_check_coords_1(self):
        dir='1'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5)
        self.assertEqual(len(result),1)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],67)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercee']['segtype'],'lipid')
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['segname'],'S2')
        self.assertEqual(p['piercer']['segtype'],'lipid')

    @pytest.mark.slow
    def test_ring_check_coords_2(self):
        dir='2'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5)
        self.assertEqual(len(result),2)
        p=result[1]
        self.assertEqual(p['piercee']['resid'],45)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercer']['resid'],595)
        self.assertEqual(p['piercer']['segname'],'S2')
        p=result[0]
        self.assertEqual(p['piercee']['resid'],67)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['segname'],'S2')

    @pytest.mark.slow
    def test_ring_check_coords_3(self):
        dir='3'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5)
        self.assertEqual(len(result),3)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],67)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['segname'],'S2')
        p=result[1]
        self.assertEqual(p['piercee']['resid'],45)
        self.assertEqual(p['piercee']['segname'],'S2')
        self.assertEqual(p['piercer']['resid'],595)
        self.assertEqual(p['piercer']['segname'],'S2')
        p=result[2]
        self.assertEqual(p['piercee']['resid'],243)
        self.assertEqual(p['piercee']['segname'],'S1')
        self.assertEqual(p['piercer']['resid'],529)
        self.assertEqual(p['piercer']['segname'],'S1')

    def test_ring_check_coords_5(self):
        # synthetic: BGLC pyranose ring (glycan) in XY plane, CHL1 bond (lipid) along Z piercing through center
        dir='5'
        pdb=os.path.join(dir,'test.pdb')
        psf=os.path.join(dir,'test.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5,segtypes=['lipid','glycan'])
        self.assertEqual(len(result),1)
        p=result[0]
        self.assertEqual(p['piercee']['segname'],'AG01')
        self.assertEqual(p['piercee']['resid'],1)
        self.assertEqual(p['piercee']['segtype'],'glycan')
        self.assertEqual(p['piercer']['segname'],'LIP1')
        self.assertEqual(p['piercer']['resid'],2)
        self.assertEqual(p['piercer']['segtype'],'lipid')

    def test_ring_check_resname_and_only_piercees(self):
        # new fields/filters used by the protein side-chain auto-rotate path
        dir='5'
        pdb=os.path.join(dir,'test.pdb')
        psf=os.path.join(dir,'test.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5,segtypes=['lipid','glycan'])
        self.assertEqual(len(result),1)
        # resname is now reported and actually resolved (needed to name the fix, e.g.
        # 'HIS G-330'); it must not fall back to '?' (a resid int/str key mismatch bug)
        self.assertNotEqual(result[0]['piercee'].get('resname'), '?')
        self.assertNotEqual(result[0]['piercer'].get('resname'), '?')
        self.assertTrue(result[0]['piercee']['resname'])
        # only_piercees restricts which rings are checked (targeted re-check after a rotation)
        hit=ring_check(psf,pdb,xsc,cutoff=3.5,segtypes=['lipid','glycan'],only_piercees=[('AG01',1)])
        self.assertEqual(len(hit),1)
        miss=ring_check(psf,pdb,xsc,cutoff=3.5,segtypes=['lipid','glycan'],only_piercees=[('ZZZZ',999)])
        self.assertEqual(len(miss),0)

    @pytest.mark.slow
    def test_ring_check_coords_4(self):
        # checks when molecules are in different periodic images
        dir='4'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc,cutoff=3.5)
        self.assertEqual(len(result),1)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],907)
        self.assertEqual(p['piercee']['segname'],'S1')
        self.assertEqual(p['piercer']['resid'],354)
        self.assertEqual(p['piercer']['segname'],'S2')