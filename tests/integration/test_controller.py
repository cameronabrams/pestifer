"""

.. module:: test_controller
   :synopsis: tests pestifer.Controller
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
from pestifer.controller import Controller
from pestifer.resourcemanager import ResourceManager

class TestController(unittest.TestCase):
    def setUp(self):
        self.userinputs='example.yaml'
    def test_controller_init(self):
        c=Controller(self.userinputs)
        expected_namd2='/alt/path/to/namd2'
        self.assertEqual(expected_namd2,ResourceManager.namd2)
