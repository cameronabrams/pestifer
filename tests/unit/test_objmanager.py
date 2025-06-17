# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
from pestifer.objmanager import ObjManager

class TestObjManager(unittest.TestCase):
    def test_obj_manager(self):
        O=ObjManager()
        self.assertEqual(len(O.obj_classes),16)