import unittest
from pestifer.config import Config
from pestifer.controller import Controller
from pestifer.resourcemanager import ResourceManager

class TestController(unittest.TestCase):
    def test_controller_base(self):
        RM=ResourceManager()
        configfile=RM.get_example_yaml_by_index(1)
        config=Config(configfile)
        C=Controller(config,index=1)
        self.assertEqual(C.config['user']['title'],'BPTI')
        self.assertEqual(len(C.tasks),12)
        self.assertEqual(C.index,1)
        # configfile=RM.get_example_yaml_by_index(7)
        # C=Controller(configfile)
        # self.assertEqual(C.config['user']['title'],'HIV-1 Env Trimer 8fad, drug molecule removed')
        # self.assertEqual(len(C.tasks),15)
