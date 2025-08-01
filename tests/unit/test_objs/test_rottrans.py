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

class TestRotTransList(unittest.TestCase):

    def test_rottrans_list_creation(self):
        rottrans_list = RotTransList()
        self.assertIsInstance(rottrans_list, RotTransList)
