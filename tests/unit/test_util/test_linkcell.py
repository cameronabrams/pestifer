from pestifer.util.linkcell import Linkcell
from pestifer.util.util import cell_from_xsc
import unittest
import numpy as np

class TestLinkcell(unittest.TestCase):
    def test_linkcell_create(self):
        box,orig=cell_from_xsc('test.xsc')
        sidelengths=np.diagonal(box)
        ll=orig-0.5*sidelengths
        ur=orig+0.5*sidelengths
        LC=Linkcell(np.array([ll,ur]),10.0)
        self.assertEqual(LC.cells_per_dim[0],8)
        self.assertEqual(LC.cells_per_dim[1],8)
        self.assertEqual(LC.cells_per_dim[2],11)

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
    
    def test_linkcell_indexing(self):
        box=np.array([[100,0,0],[0,100,0],[0,0,100]],dtype=float)
        orig=np.array([0,0,0],dtype=float)
        sidelengths=np.diagonal(box)
        ll=orig-0.5*sidelengths
        ur=orig+0.5*sidelengths
        LC=Linkcell(np.array([ll,ur]),10.0)
        for p in np.random.random(3000).reshape(1000,3)*(100)-52:
            c=LC.cellndx_of_point(p)
            lc=LC.ldx_of_cellndx(c)
            cc=LC.cellndx_of_ldx(lc)
            # logger.debug(f'{c} {lc} {cc}')
            assert all(c==cc)
    # def test_populate(self):
    #     box,orig=cell_from_xsc('test.xsc')
    #     sidelengths=np.diagonal(box)
    #     ll=orig-0.5*sidelengths
    #     ur=orig+0.5*sidelengths
    #     LC=Linkcell(np.array([ll,ur]),10.0,atidxlabel=None)
    #     coorddf=coorddf_from_pdb('test')
    #     self.assertFalse(coorddf.isnull().values.any())
    #     LC.populate(coorddf,ncpu=os.cpu_count())
    #     self.assertEqual(np.round(LC.avg_cell_pop,0),109.0)
    #     self.assertTrue('linkcell_idx' in coorddf)

