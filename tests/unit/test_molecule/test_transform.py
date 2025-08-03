import unittest
from pestifer.molecule.transform import Transform, TransformList
from pidibble.pdbparse import PDBRecord, PDBRecordList
import numpy as np
from dataclasses import dataclass

class TestTransform(unittest.TestCase):

    def test_transform_identity(self):
        bi = Transform.identity()
        I = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=float)
        self.assertTrue(np.array_equal(I, bi.tmat))

    def test_transform_from_matvec(self):
        R = np.array([[1, 1, 1], [0, 1, 0], [0, 0, 1]], dtype=float)
        RT = np.array([[1, 1, 1, 5], [0, 1, 0, 4], [0, 0, 1, 3], [0, 0, 0, 1]], dtype=float)
        T = np.array([5, 4, 3], dtype=float)
        bi = Transform(R, T)
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
        bi = Transform(barec)
        self.assertIsNotNone(bi.tmat)

class TestTransformList(unittest.TestCase):

    def test_transform_list_initialization(self):
        tl = TransformList()
        self.assertIsInstance(tl, TransformList)
        self.assertEqual(len(tl), 0)

    def test_transform_list_from_transforms(self):
        transforms = [Transform.identity(), Transform.identity()]
        tl = TransformList(transforms)
        self.assertEqual(len(tl), 2)
        self.assertIsInstance(tl[0], Transform)

    def test_transform_list_from_pdb_records(self):
        @dataclass
        class test_row:
            m1: float = 0.0
            m2: float = 0.0
            m3: float = 0.0
            t: float = 1.0

        barec1 = PDBRecord(input_dict={'row': [test_row(m1=1.0), test_row(m2=1.0), test_row(m3=1.0)], 'coordinate': [1, 2, 3]})
        barec2 = PDBRecord(input_dict={'row': [test_row(m1=1.0), test_row(m2=1.0), test_row(m3=1.0)], 'coordinate': [1, 2, 3]})
        tl = TransformList(PDBRecordList([barec1, barec2]))
        self.assertEqual(len(tl), 2)
        self.assertIsInstance(tl[0], Transform)


