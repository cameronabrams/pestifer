import unittest
from pestifer.controller import Controller
from pestifer.mods import *
from pestifer.config import ResourceManager
class TestController(unittest.TestCase):
    def test_controller_base(self):
        rm=ResourceManager()
        configfile=rm['examples']+'/01-bpti.yaml'
        C=Controller(configfile)
        self.assertEqual(C.config['user']['title'],'BPTI')
        self.assertEqual(len(C.tasks),4)
        configfile=rm['examples']+'/06-hiv-env-8fad.yaml'
        C=Controller(configfile)
        self.assertEqual(C.config['user']['title'],'HIV-1 Env Trimer 8fad, drug molecule removed')
        self.assertEqual(len(C.tasks),6)
