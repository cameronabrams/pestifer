from pestifer.linkcell import Linkcell
import unittest
from pestifer.tasks import *
import os
from collections import UserList
import numpy as np


class TestLinkcell(unittest.TestCase):
    def test_linkcell(self):
        box,orig=cell_from_xsc('test.xsc')
        LC=Linkcell(box,10.0,origin=orig)
        self.assertEqual(LC.ncells[0],8)
        self.assertEqual(LC.ncells[1],8)
        self.assertEqual(LC.ncells[2],11)

    def test_wrap(self):
        box=np.array([[100,0,0],[0,100,0],[0,0,100]],dtype=float)
        orig=np.array([0,0,0],dtype=float)
        LC=Linkcell(box,10.0,origin=orig)
        p=np.array([-2,-2,-2])
        pp,bl=LC.wrap_point(p)
        self.assertEqual(pp[0],98)
        self.assertEqual(pp[1],98)
        self.assertEqual(pp[2],98)
        self.assertEqual(bl[0],1)
    
    def test_populate(self):
        from pidibble.pdbparse import PDBParser
        box,orig=cell_from_xsc('test.xsc')
        LC=Linkcell(box,10.0,origin=orig)
        p=PDBParser(PDBcode='test').parse()
        atlist=p.parsed['ATOM']
        serial=[x.serial for x in atlist]
        x=[x.x for x in atlist]
        y=[x.y for x in atlist]
        z=[x.z for x in atlist]
        coorddf=pd.DataFrame({'serial':serial,'x':x,'y':y,'z':z})
        LC.populate(coorddf,ncpu=os.cpu_count())
        self.assertEqual(np.round(LC.avg_cell_pop,0),109.0)
        self.assertTrue('linkcell_idx' in coorddf)

