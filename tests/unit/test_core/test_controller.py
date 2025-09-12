import logging
import os
import unittest
from pestifer.tasks.basetask import BaseTask
from pestifer.core.config import Config
from pestifer.core.controller import Controller
from pestifer.core.resourcemanager import ResourceManager

logger = logging.getLogger(__name__)    

class TestController(unittest.TestCase):

    def test_controller_base(self):
        RM = ResourceManager()
        EM = RM.example_manager
        # read in config file directly from the example in the package
        e1 = EM.checkout_example(1)
        config = Config(userfile=e1.scriptname).configure_new()
        C = Controller().configure(config, userspecs={'title': 'Bovine Pancreatic Trypsin Inhibitor (BPTI)'}, index=1)

        self.assertEqual(C.config['user']['title'],'Bovine Pancreatic Trypsin Inhibitor (BPTI)')
        self.assertEqual(len(C.tasks), 14)
        self.assertEqual(C.index, 1)
        task1 = C.tasks[0]
        self.assertEqual(task1.taskname, 'fetch')
        self.assertTrue('sourceID' in task1.specs)
        self.assertEqual(task1.specs['sourceID'], '6pti')
        task2 = C.tasks[1]
        self.assertEqual(task2.taskname, 'psfgen')
        task3 = C.tasks[2]
        self.assertEqual(task3.taskname, 'validate')
        self.assertTrue(all([x.is_provisioned for x in C.tasks]))

    def test_controller_spawn_subcontroller(self):
        config=Config().configure_new()
        self.assertEqual(len(config['user']['tasks']),0)
        config['user']['tasks'] = [{'fetch': {'sourceID': '6pti'}}]
        C = Controller().configure(config, userspecs={'title': 'Test Controller'}, index=0)
        self.assertTrue(C.tasks[0].is_provisioned)
        subC = Controller.spawn_subcontroller(C)
        self.assertEqual(subC.index, 1)
        C.tasks[0].subcontroller = subC
        subcontroller = C.tasks[0].subcontroller
        self.assertIsInstance(subcontroller, Controller)
        self.assertEqual(subcontroller.index, 1)
        self.assertNotEqual(subcontroller.config, C.config)
        self.assertEqual(len(subcontroller.tasks), 0)
        namd_global_config = C.tasks[0].provisions.get('namd_global_config', {})
        self.assertTrue('namd-version' in namd_global_config)
        # Now reconfigure the subcontroller with new tasks
        new_tasks = [{'fetch': {}}, {'psfgen': {}}, {'md': {}}]
        subcontroller.reconfigure_tasks(new_tasks)
        self.assertEqual(len(subcontroller.tasks), 3)
        self.assertEqual(subcontroller.tasks[0].taskname, 'fetch')
        self.assertEqual(subcontroller.tasks[1].taskname, 'psfgen')
        self.assertEqual(subcontroller.tasks[2].taskname, 'md')
        self.assertTrue(all([x.is_provisioned for x in subcontroller.tasks]))