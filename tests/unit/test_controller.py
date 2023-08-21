import unittest
from pestifer.controller import Controller
from pestifer.mods import *

class TestController(unittest.TestCase):
    def test_controller_base(self):
        C=Controller('user_config.yaml')
        self.assertEqual(C.config.defs['title'],'Test User Config')
        self.assertEqual(len(C.tasks),2)

    def test_controller_do(self):
        C=Controller('user_config.yaml')
        C.do_tasks()
        # self.assertTrue(os.path.exists('step1-layloops.pdb'))
