import os
import unittest
from pestifer.tasks import task_classes
from pestifer.tasks.basetask import BaseTask
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
        C=Controller().configure(config, userspecs={'title': 'Bovine Pancreatic Trypsin Inhibitor (BPTI)'}, index=1)

        self.assertEqual(C.config['user']['title'],'Bovine Pancreatic Trypsin Inhibitor (BPTI)')
        self.assertEqual(len(C.tasks),13)
        self.assertEqual(C.index,1)
        task1 = C.tasks[0]
        self.assertEqual(task1.taskname, 'fetch')
        self.assertTrue('sourceID' in task1.specs)
        self.assertEqual(task1.specs['sourceID'], '6pti')
        task2 = C.tasks[1]
        self.assertEqual(task2.taskname, 'psfgen')
        task3 = C.tasks[2]
        self.assertEqual(task3.taskname, 'md')
        self.assertTrue(all([x.is_provisioned for x in C.tasks]))

    def test_controller_spawn_subcontroller(self):
        RM=ResourceManager()
        config=Config()
        class FakeTask(BaseTask):
            def __init__(self, specs=None, provisions=None):
                super().__init__(specs, provisions)
            def do(self):
                return {'status': 'done'}

        task_classes['fake_task'] = FakeTask

        config['user']['tasks'] = [{'fake_task': {}}]
        C = Controller().configure(config, userspecs={'title': 'Test Controller'}, index=0)
        self.assertTrue(C.tasks[0].is_provisioned)
        subcontroller = C.tasks[0].subcontroller
        self.assertIsInstance(subcontroller, Controller)
        self.assertEqual(subcontroller.index, 1)
        self.assertEqual(subcontroller.config, C.config)
        self.assertEqual(subcontroller.tasks, None)
        namd_global_config = C.tasks[0].provisions.get('namd_global_config', {})
        self.assertTrue('namd-version' in namd_global_config)
        # Now reconfigure the subcontroller with new tasks
        new_tasks = [{'fake_task': {'taskname': 'subtask1'}}, {'fake_task': {'taskname': 'subtask2'}}]
        subcontroller.reconfigure_tasks(new_tasks)
        self.assertEqual(len(subcontroller.tasks), 2)
        self.assertEqual(subcontroller.tasks[0].taskname, 'subtask1')
        self.assertEqual(subcontroller.tasks[1].taskname, 'subtask2')
        self.assertTrue(all([x.is_provisioned for x in subcontroller.tasks]))