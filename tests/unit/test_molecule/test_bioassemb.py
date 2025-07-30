import unittest
from pestifer.molecule.bioassemb import BioAssemb, Transform, TransformList, BioAssembList
from mmcif.api.PdbxContainers import DataContainer
from pestifer.molecule.asymmetricunit import AsymmetricUnit
from pidibble.pdbparse import PDBRecord, PDBRecordDict
import numpy as np
from dataclasses import dataclass

class TestTransform(unittest.TestCase):

    def test_transform_empty(self):
        bi = Transform()
        self.assertEqual(bi.tmat, None)
        self.assertEqual(bi.index, 0)
        self.assertEqual(bi.applies_chainIDs, None)
        self.assertEqual(bi.chainIDmap, None)
        self.assertEqual(bi.segname_by_type_map, None)

    def test_transform_identity(self):
        bi = Transform.identity()
        I = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=float)
        self.assertTrue(np.array_equal(I, bi.tmat))

    def test_transform_from_matvec(self):
        R = np.array([[1, 1, 1], [0, 1, 0], [0, 0, 1]], dtype=float)
        RT = np.array([[1, 1, 1, 5], [0, 1, 0, 4], [0, 0, 1, 3], [0, 0, 0, 1]], dtype=float)
        T = np.array([5, 4, 3], dtype=float)
        bi = Transform.from_matvec(R, T)
        self.assertTrue(np.array_equal(RT, bi.tmat))

    def test_transform_from_pdbrecord(self):
        @dataclass
        class test_row:
            m1: float = 0.0
            m2: float = 0.0
            m3: float = 0.0
            t: float = 1.0
        input_dict = {'row': [test_row(m1=1.0), test_row(m2=1.0), test_row(m3=1.0)],
                      'coordinate': [1, 2, 3]}
        barec = PDBRecord(input_dict=input_dict)
        bi = Transform.from_pdb_record(barec)
        self.assertIsNotNone(bi.tmat)



# class TestBiomT(unittest.TestCase):
#     def test_identity(self):
#         bi=Transform('identity')
#         I=np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0],[0, 0, 0, 1]],dtype=float)
#         self.assertTrue(np.array_equal(I,bi.tmat))
#     def test_actual(self):
#         R=np.array([[1, 1, 1],[0, 1, 0],[0, 0, 1]],dtype=float)
#         RT=np.array([[1, 1, 1, 5],[0, 1, 0, 4],[0, 0, 1, 3],[0, 0, 0, 1]],dtype=float)
#         T=np.array([5, 4, 3],dtype=float)
#         bi=Transform(R,T,[],0)
#         self.assertTrue(np.array_equal(RT,bi.tmat))
#     def test_au(self):
#         au=AsymmetricUnit() # empty
#         ba=BioAssemb(au)
#         bi=ba.transforms[0]
#         I=np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0],[0, 0, 0, 1]],dtype=float)
#         self.assertTrue(np.array_equal(I,bi.tmat))

