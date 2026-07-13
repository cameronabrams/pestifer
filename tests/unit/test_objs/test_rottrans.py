import unittest

from pestifer.objs.rottrans import RotTrans, RotTransList

class TestRotTrans(unittest.TestCase):

    def test_rottrans_creation(self):
        rottrans = RotTrans(movetype='TRANS', x=1.0, y=2.0, z=3.0)
        self.assertEqual(rottrans.movetype, 'TRANS')
        self.assertEqual(rottrans.x, 1.0)
        self.assertEqual(rottrans.y, 2.0)
        self.assertEqual(rottrans.z, 3.0)
        self.assertIsInstance(rottrans, RotTrans)

    def test_rottrans_from_shortcode(self):
        rottrans = RotTrans('TRANS,1.0,2.0,3.0')
        self.assertIsInstance(rottrans, RotTrans)
        self.assertEqual(rottrans.movetype, 'TRANS')
        self.assertEqual(rottrans.x, 1.0)
        self.assertEqual(rottrans.y, 2.0)
        self.assertEqual(rottrans.z, 3.0)
        rottrans = RotTrans('ROT,x,90.0')
        self.assertIsInstance(rottrans, RotTrans)
        self.assertEqual(rottrans.movetype, 'ROT')
        self.assertEqual(rottrans.axis, 'x')
        self.assertEqual(rottrans.angle, 90.0)

    def test_axisangle_creation(self):
        # AXISANGLE is the rigid-body rotation formerly offered as the ANGLEIJK irotation:
        # rotate a fragment about the axis normal to a three-atom angle, pivoting at selJ.
        atoms = ['segname A and resid 12 and name N',
                 'segname B and resid 15 and name C',
                 'segname B and resid 18 and name O']
        rt = RotTrans(movetype='AXISANGLE', axis_atoms=atoms, angle=180.0, sel='segname B')
        self.assertEqual(rt.movetype, 'AXISANGLE')
        self.assertEqual(rt.axis_atoms, atoms)
        self.assertEqual(rt.angle, 180.0)
        self.assertEqual(rt.sel, 'segname B')

    def test_axisangle_requires_three_atomselections(self):
        with self.assertRaises(ValueError):
            RotTrans(movetype='AXISANGLE', axis_atoms=['a', 'b'], angle=90.0)  # only two
        with self.assertRaises(ValueError):
            RotTrans(movetype='AXISANGLE', angle=90.0)  # missing axis_atoms

class TestRotTransList(unittest.TestCase):

    def test_rottrans_list_creation(self):
        rottrans_list = RotTransList()
        self.assertIsInstance(rottrans_list, RotTransList)
