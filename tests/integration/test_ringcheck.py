import unittest
from pestifer.linkcell import Linkcell
from pestifer.ring import RingList, ring_check
from pestifer.coord import coorddf_from_pdb,mic_shift
from pestifer.util import cell_from_xsc
from pestifer.psf import PSFContents,PSFBondList
import os
import logging
from itertools import compress
logger=logging.getLogger(__name__)
import numpy as np

class TestRingCheck(unittest.TestCase):
    def test_ring_check_coords_1(self):
        dir='1'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc)
        self.assertEqual(len(result),1)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],67)
        self.assertEqual(p['piercee']['chain'],'S2')
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['chain'],'S2')

    def test_ring_check_coords_2(self):
        dir='2'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc)
        self.assertEqual(len(result),1)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],67)
        self.assertEqual(p['piercee']['chain'],'S2')
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['chain'],'S2')

    def test_ring_check_coords_3(self):
        dir='3'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc)
        self.assertEqual(len(result),2)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],67)
        self.assertEqual(p['piercee']['chain'],'S2')
        self.assertEqual(p['piercer']['resid'],645)
        self.assertEqual(p['piercer']['chain'],'S2')
        p=result[1]
        self.assertEqual(p['piercee']['resid'],243)
        self.assertEqual(p['piercee']['chain'],'S1')
        self.assertEqual(p['piercer']['resid'],529)
        self.assertEqual(p['piercer']['chain'],'S1')

    def test_ring_check_coords_4(self):
        # checks when molecules are in different periodic images
        dir='4'
        pdb=os.path.join(dir,'S2.pdb')
        psf=os.path.join(dir,'S2.psf')
        xsc=os.path.join(dir,'test.xsc')
        result=ring_check(psf,pdb,xsc)
        self.assertEqual(len(result),1)
        p=result[0]
        self.assertEqual(p['piercee']['resid'],907)
        self.assertEqual(p['piercee']['chain'],'S1')
        self.assertEqual(p['piercer']['resid'],354)
        self.assertEqual(p['piercer']['chain'],'S2')