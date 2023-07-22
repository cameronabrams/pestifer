"""

.. module:: test_controller
   :synopsis: tests pestifer.Controller
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import unittest
import platform
import os
from pestifer.controller import Controller

class TestController(unittest.TestCase):
    def setUp(self):
        self.userinputs='example.yaml'
        plat=platform.system()
        if plat=='Linux':
            self.user_home=os.environ['HOME']
        elif plat=='Windows':
            self.user_home=os.environ['HOMEPATH']
        elif plat=='Darwin':
            self.user_home=os.environ['HOME']
        else:
            raise(Exception,f'platform {plat} not recognized')
    def test_controller_init(self):
        c=Controller(self.userinputs)
        expected_namd2='/alt/path/to/namd2'
        self.assertEqual(expected_namd2,c.resman.namd2)