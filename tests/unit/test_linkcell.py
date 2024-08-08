from pestifer.linkcell import Linkcell
from pestifer.util import cell_from_xsc
import unittest
from pestifer.tasks import *
import os
from collections import UserList
import numpy as np


class TestLinkcell(unittest.TestCase):
    def test_linkcell_create(self):
        box,orig=cell_from_xsc('test.xsc')
        sidelengths=np.diagonal(box)
        ll=orig-0.5*sidelengths
        ur=orig+0.5*sidelengths
        LC=Linkcell(np.array([ll,ur]),10.0)
        self.assertEqual(LC.ncells[0],8)
        self.assertEqual(LC.ncells[1],8)
        self.assertEqual(LC.ncells[2],11)

    def test_wrap(self):
        box=np.array([[100,0,0],[0,100,0],[0,0,100]],dtype=float)
        orig=np.array([0,0,0],dtype=float)
        sidelengths=np.diagonal(box)
        ll=orig-0.5*sidelengths
        ur=orig+0.5*sidelengths
        LC=Linkcell(np.array([ll,ur]),10.0)
        
        p=np.array([-52,-52,-52])
        pp,bl=LC.wrap_point(p)
        self.assertEqual(pp[0],48)
        self.assertEqual(pp[1],48)
        self.assertEqual(pp[2],48)
        self.assertEqual(bl[0],1)
    
    def test_populate(self):
        from pidibble.pdbparse import PDBParser
        box,orig=cell_from_xsc('test.xsc')
        sidelengths=np.diagonal(box)
        ll=orig-0.5*sidelengths
        ur=orig+0.5*sidelengths
        LC=Linkcell(np.array([ll,ur]),10.0)
        p=PDBParser(PDBcode='test').parse()
        atlist=p.parsed['ATOM']
        serial=[x.serial for x in atlist]
        x=[x.x for x in atlist]
        y=[x.y for x in atlist]
        z=[x.z for x in atlist]
        if not all([_==_ for _ in x]):
            for i,_ in enumerate(x):
                if not _==_:
                    logger.debug(f'{atlist[i]}')
        self.assertTrue(all([_==_ for _ in y]))
        self.assertTrue(all([_==_ for _ in z]))
        coorddf=pd.DataFrame({'serial':serial,'x':x,'y':y,'z':z})
        self.assertFalse(coorddf.isnull().values.any())
        LC.populate(coorddf,ncpu=os.cpu_count())
        self.assertEqual(np.round(LC.avg_cell_pop,0),109.0)
        self.assertTrue('linkcell_idx' in coorddf)

