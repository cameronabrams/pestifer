import unittest
import shutil
import os
from pestifer.core.examplemanager import ExampleManager

class TestExampleManager(unittest.TestCase):

    def setUp(self):
        if os.path.isdir('./package'):
            shutil.rmtree('./package', ignore_errors=True)
        if os.path.isdir('./userspace'):
            shutil.rmtree('./userspace', ignore_errors=True)
        shutil.copytree('../fixtures/test_example_manager_inputs/package', './package/', dirs_exist_ok=True)
        shutil.copytree('../fixtures/test_example_manager_inputs/userspace', './userspace/', dirs_exist_ok=True)
        self.manager = ExampleManager('package/package/resources/examples',docs_source_path='package/docs/source')

    def tearDown(self):
        shutil.rmtree('./package', ignore_errors=True)
        shutil.rmtree('./userspace', ignore_errors=True)

    def test_example_manager_init(self):
        self.assertIsNotNone(self.manager)
        self.assertEqual(len(self.manager.examples_list), 2)
        self.assertEqual(self.manager.examples_list[0].name, 'exA')
        self.assertEqual(self.manager.examples_list[1].name, 'exB')
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'info.yaml')))
        self.assertTrue(os.path.isdir(self.manager.path))

    def test_example_manager_new_example(self):
        os.chdir('userspace')
        example = self.manager.new_example('exX.yaml')
        self.assertTrue(example.description.startswith('exX: Use Pestifer to create a simple userspace application'))
        self.assertEqual(example.pdbID, '1a8l')
        self.assertEqual(example.name, 'exX')
        self.assertEqual(example.index, 0)  # Index should be set to 0
        os.chdir('..')

    def test_example_manager_add_example(self):
        os.chdir('userspace')
        newidx=self.manager.add_example('exX.yaml')
        os.chdir('..')
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'exX.yaml')))
        self.assertEqual(newidx, 3)

    def test_example_manager_insert_example(self):
        os.chdir('userspace')
        reqidx=self.manager.insert_example(2,'exX')
        os.chdir('..')
        self.assertEqual(reqidx, 2)
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'exX.yaml')))
        self.assertEqual(self.manager.examples_list[1].name, 'exX')
    
    def test_example_manager_delete_example(self):
        self.manager.delete_example(2)
        self.assertFalse(os.path.isfile(os.path.join(self.manager.path, 'exB.yaml')))

    def test_example_manager_update_example(self):
        os.chdir('userspace')
        self.manager.update_example(1, 'exA_updated.yaml')
        os.chdir('..')
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'exA_updated.yaml')))
        self.assertEqual(self.manager.examples_list[0].name, 'exA_updated')