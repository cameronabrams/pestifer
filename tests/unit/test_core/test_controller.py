import os
import unittest
from pestifer.core.config import Config
from pestifer.core.controller import Controller
from pestifer.core.resourcemanager import ResourceManager

class TestController(unittest.TestCase):
    def test_controller_base(self):
        RM=ResourceManager()
        EM=RM.example_manager
        # read in config file directly from the example in the package
        configfile=os.path.join(EM.path,EM.examples_list[0].name,EM.examples_list[0].name+'.yaml')
        config=Config(configfile)
        C=Controller(config,index=1)
        self.assertEqual(C.config['user']['title'],'Bovine Pancreatic Trypsin Inhibitor (BPTI)')
        self.assertEqual(len(C.tasks),12)
        self.assertEqual(C.index,1)
