# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import unittest
from pestifer.core.objmanager import ObjManager

class TestObjManager(unittest.TestCase):
    def test_obj_manager(self):
        O=ObjManager()
        self.assertEqual(len(O._obj_classes),16)