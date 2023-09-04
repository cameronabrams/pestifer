import unittest
from pestifer.controller import Controller
from pestifer.mods import *
from pestifer.config import ResourceManager
class TestController(unittest.TestCase):
    def test_controller_base(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/example1.yaml'
        C=Controller(configfile)
        print(C.config['base'])
        self.assertEqual(C.config['user']['title'],'BPTI')
        self.assertEqual(len(C.tasks),4)
        configfile=rm['examples']+'/example6.yaml'
        C=Controller(configfile)
        print(C.config['base'])
        self.assertEqual(C.config['user']['title'],'HIV-1 Env Trimer 8fad')
        self.assertEqual(len(C.tasks),6)
