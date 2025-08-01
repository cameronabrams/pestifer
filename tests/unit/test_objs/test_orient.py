import unittest
from pestifer.objs.orient import Orient, OrientList

class TestOrient(unittest.TestCase):
    def test_orient_creation(self):
        orient = Orient(
            axis="x",
            refatom="CA"
        )
        self.assertIsInstance(orient, Orient)
        self.assertEqual(orient.axis, "x")
        self.assertEqual(orient.refatom, "CA")
        self.assertEqual(orient._yaml_header, 'orient')
        self.assertEqual(orient._objcat, 'coord')
        self.assertEqual(repr(orient), "Orient(axis='x', refatom='CA')")

    def test_orient_copy(self):
        orient1 = Orient(
            axis="y",
            refatom="CB"
        )
        orient2 = orient1.copy()
        self.assertTrue(orient1 is not orient2)
        self.assertEqual(orient1, orient2)
        self.assertIsInstance(orient2, Orient)
        self.assertEqual(orient2.axis, "y")
        self.assertEqual(orient2.refatom, "CB")