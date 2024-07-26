from pestifer.ring import RingList
from pestifer.util import *
import unittest
from pestifer.tasks import *
from pestifer.psf import PSFContents
import os
from collections import UserList
import numpy as np


class TestRings(unittest.TestCase):
    def test_ringlist(self):
        self.psf=PSFContents('test.psf',parse_topology=True)
        self.G=self.psf.bonds.to_graph()
        self.Rings=RingList(self.G)
        logger.debug(f'{len(self.Rings)} rings detected.')
        nrings=len(self.Rings)
        self.assertEqual(nrings,42)
        nrings_by_size={}
        for r in self.Rings:
            s=len(r)
            if not s in nrings_by_size:
                nrings_by_size[s]=0
            nrings_by_size[s]+=1
        self.assertTrue(6 in nrings_by_size)
        self.assertTrue(5 in nrings_by_size)
        self.assertFalse(4 in nrings_by_size)
        self.assertFalse(7 in nrings_by_size)

    def test_ringlist_coords(self):
        self.box,self.orig=cell_from_xsc('test.xsc')
        self.coorddf=coorddf_from_pdb('test.pdb')
        self.topol=PSFContents('test.psf',parse_topology=True)
        self.Rings=RingList(self.topol.G)
        self.Rings.injest_coordinates(self.coorddf,idx_key='serial',pos_key=['x','y','z'],box=self.box)
        self.assertTrue(all([x.same_image for x in self.Rings]))