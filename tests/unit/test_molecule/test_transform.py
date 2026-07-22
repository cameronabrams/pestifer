import unittest
from pestifer.molecule.transform import Transform, TransformList
from pestifer.molecule.atom import Atom, AtomList
from pestifer.objs.resid import ResID
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

class TestTransformApply(unittest.TestCase):

    def test_apply_array_matches_matrix_multiply(self):
        R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        tf = Transform(R, np.array([1.0, 2.0, 3.0]))
        pts = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        out = tf.apply(pts)
        self.assertTrue(np.allclose(out, [[1.0, 3.0, 3.0], [0.0, 2.0, 3.0]]))

    def test_apply_array_does_not_mutate_input(self):
        tf = Transform(np.identity(3), np.array([5.0, 0.0, 0.0]))
        pts = np.array([[1.0, 1.0, 1.0]])
        orig = pts.copy()
        tf.apply(pts)
        self.assertTrue(np.array_equal(pts, orig))

    def test_identity_apply_is_noop(self):
        tf = Transform.identity()
        pts = np.array([[1.0, 2.0, 3.0]])
        self.assertTrue(np.allclose(tf.apply(pts), pts))

    def test_apply_mutates_atomlist_in_place(self):
        a1 = Atom(serial=1, name='C', altloc=' ', resname='ALA', chainID='A', resid=ResID(1),
                  x=1.0, y=0.0, z=0.0, occ=1.0, beta=0.0, elem='C', charge=' ')
        a2 = Atom(serial=2, name='O', altloc=' ', resname='ALA', chainID='A', resid=ResID(1),
                  x=0.0, y=1.0, z=0.0, occ=1.0, beta=0.0, elem='O', charge=' ')
        al = AtomList([a1, a2])
        tf = Transform(np.identity(3), np.array([10.0, 20.0, 30.0]))
        ret = tf.apply(al)
        self.assertIs(ret, al)
        self.assertEqual((a1.x, a1.y, a1.z), (11.0, 20.0, 30.0))
        self.assertEqual((a2.x, a2.y, a2.z), (10.0, 21.0, 30.0))


class TestTransformSuperpose(unittest.TestCase):

    def test_superpose_recovers_transform_and_applies(self):
        rng = np.random.default_rng(0)
        P = rng.normal(size=(10, 3))
        R = np.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
        t = np.array([2.0, -1.0, 4.0])
        Q = P @ R.T + t
        tf, rmsd = Transform.superpose(P, Q)
        self.assertAlmostEqual(rmsd, 0.0, places=9)
        self.assertTrue(np.allclose(tf.apply(P), Q, atol=1e-9))

    def test_superpose_accepts_atomlists(self):
        def al(coords):
            atoms = [Atom(serial=i + 1, name='C', altloc=' ', resname='ALA', chainID='A',
                          resid=ResID(1), x=x, y=y, z=z, occ=1.0, beta=0.0, elem='C', charge=' ')
                     for i, (x, y, z) in enumerate(coords)]
            return AtomList(atoms)
        P = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        Q = P + np.array([5.0, 6.0, 7.0])
        tf, rmsd = Transform.superpose(al(P), al(Q))
        self.assertAlmostEqual(rmsd, 0.0, places=9)
        self.assertTrue(np.allclose(tf.apply(P), Q, atol=1e-9))


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


