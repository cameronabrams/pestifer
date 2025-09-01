from pestifer.tasks.cleave import CleaveTask
from pestifer.core.controller import Controller
from pestifer.core.config import Config
import os
import unittest
class TestCleaveTask(unittest.TestCase):
    def setUp(self):
        self.controller = Controller().configure(Config().configure_new(), terminate=False)
        self.scripters = self.controller.config.scripters # shortcut

    def test_cleave(self):
       tasklist = [
           {'continuation': {'psf':'in.psf','pdb':'in.pdb'}},
           {'cleave': { 
               'sites': [
                   'A:685-686',
                   'B:685-686',
                   'C:685-686'
               ]
           }}
       ]
       self.controller.reconfigure_tasks(tasklist)
       result = self.controller.do_tasks()
       self.assertTrue(os.path.exists('00-01-00_cleave-build.pdb'))
