import unittest
from pestifer.ring import ring_check
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
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['segname'],'S2')

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