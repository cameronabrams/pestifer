import unittest
from pestifer.psfutil.psfbond import PSFBond
from pestifer.psfutil.psfangle import PSFAngle
from pestifer.psfutil.psfdihedral import PSFDihedral


class TestPSFTopoelements(unittest.TestCase):

    def test_psf_bond_eq(self):
        b1 = PSFBond([1, 2])
        b2 = PSFBond([2, 1])
        self.assertEqual(b1, b2)

    def test_psf_angle_eq(self):
        a1 = PSFAngle([1, 2, 3])
        a2 = PSFAngle([3, 2, 1])
        self.assertEqual(a1, a2)

    def test_psf_dihedral_eq(self):
        d1 = PSFDihedral([1, 2, 3, 4])
        d2 = PSFDihedral([4, 3, 2, 1])
        self.assertEqual(d1, d2)

