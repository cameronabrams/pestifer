import unittest
import shutil
import os
import yaml
from pathlib import Path
from pestifer.core.examplemanager import ExampleManager

class TestExampleManager(unittest.TestCase):

    def setUp(self):
        if os.path.isdir('./project'):
            shutil.rmtree('./project', ignore_errors=True)
        if os.path.isdir('./userspace'):
            shutil.rmtree('./userspace', ignore_errors=True)
        shutil.copytree('../fixtures/test_example_manager_inputs/project', './project/', dirs_exist_ok=True)
        shutil.copytree('../fixtures/test_example_manager_inputs/userspace', './userspace/', dirs_exist_ok=True)
        self.manager = ExampleManager(
            example_resource_folder_name = 'examples',
                          resources_path = 'project/package/resources',
                        docs_source_path = 'project/docs/source'
        )

    def tearDown(self):
        shutil.rmtree('project')
        shutil.rmtree('userspace')

    def test_example_manager_init(self):
        self.assertIsNotNone(self.manager)
        self.assertEqual(len(self.manager.examples_list), 0)
        self.assertFalse(os.path.isfile(os.path.join(self.manager.path, 'info.yaml')))
        self.assertTrue(os.path.isdir(self.manager.path))

    def _build_example_set(self):
        os.chdir('userspace')
        example = self.manager.add_example('exA.yaml', companion_files=['exA_companion1.pdb', 'exA_companion1.psf'])
        example = self.manager.add_example('exB.yaml')
        os.chdir('..')

    def test_example_manager_build_example_set(self):
        self._build_example_set()
        self.assertEqual(len(self.manager.examples_list), 2)
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'info.yaml')))
        with open(os.path.join(self.manager.path, 'info.yaml'), 'r') as f:
            info = yaml.safe_load(f)
        self.assertIn('examples', info)
        self.assertEqual(len(info['examples']), 2)
        self.assertEqual(info['examples'][0]['name'], 'exA')
        self.assertEqual(info['examples'][1]['name'], 'exB')
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'exA', 'exA.yaml')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'exA', 'exA_companion1.pdb')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'exA', 'exA_companion1.psf')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'exB', 'exB.yaml')))
        self.assertTrue(os.path.isfile(self.manager.sphinx_example_manager.examples_rst))
        self.assertTrue(os.path.isdir(self.manager.sphinx_example_manager.examples_folder_path))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.sphinx_example_manager.examples_folder_path, 'exA.rst')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.sphinx_example_manager.examples_folder_path, 'exB.rst')))

    def test_example_manager_insert_example(self):
        self._build_example_set()
        os.chdir('userspace')
        reqidx=self.manager.insert_example(2,'exX')
        os.chdir('..')
        self.assertEqual(reqidx, 2)
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, 'exX', 'exX.yaml')))
        self.assertEqual(self.manager.examples_list[1].name, 'exX')
    
    def test_example_manager_delete_example(self):
        self._build_example_set()
        self.manager.delete_example(2)
        self.assertFalse(os.path.isfile(os.path.join(self.manager.path, 'exB.yaml')))

    def test_example_manager_update_example_inplace(self):
        self._build_example_set()
        os.chdir('userspace')
        os.mkdir('exa_updated')
        os.chdir('exa_updated')
        Path('exA_companion1.pdb').touch()
        Path('exA_companion1.psf').touch()
        self.manager.update_example(1, description='Updated description', pdbID='9ggs', author_name='Mel Brooks')
        os.chdir('..')
        shutil.rmtree('exa_updated')
        os.chdir('..')
        with open(os.path.join(self.manager.path, 'info.yaml'), 'r') as f:
            info = yaml.safe_load(f)
        self.assertIn('examples', info)
        self.assertEqual(len(info['examples']), 2)
        self.assertEqual(info['examples'][0]['name'], 'exA')
        self.assertEqual(info['examples'][0]['description'], 'Updated description')
        self.assertEqual(info['examples'][0]['pdbID'], '9ggs')
        self.assertEqual(info['examples'][0]['author_name'], 'Mel Brooks')
        self.assertEqual(info['examples'][0]['index'], 1)

    def test_example_manager_update_example_reyaml(self):
        self._build_example_set()
        cwd = os.getcwd()
        os.chdir('userspace')
        os.mkdir('exa_updated')
        os.chdir('exa_updated')
        Path('exA_companion1.pdb').touch()
        Path('exA_companion1.psf').touch()
        with open('../exA.yaml', 'r') as f:
            lines= f.readlines()
        with open('exA.yaml', 'w') as f:
            for line in lines:
                if 'title:' in line:
                    f.write('title: Updated title\n')
                else:
                    f.write(line)
        self.manager.update_example(1)
        os.chdir('..')
        shutil.rmtree('exa_updated')
        with open(os.path.join(self.manager.path, 'info.yaml'), 'r') as f:
            info = yaml.safe_load(f)
        self.assertIn('examples', info)
        self.assertEqual(len(info['examples']), 2)
        self.assertEqual(info['examples'][0]['name'], 'exA')
        self.assertEqual(info['examples'][0]['description'], 'Updated title')
        self.assertEqual(info['examples'][0]['index'], 1)
        os.chdir(cwd)

    def test_example_manager_update_example_rename(self):
        cwd = os.getcwd()
        self._build_example_set()
        os.chdir('userspace')
        os.mkdir('exa_updated')
        os.chdir('exa_updated')
        Path('exA_companion1.pdb').touch()
        Path('exA_companion1.psf').touch()
        with open('../exA.yaml', 'r') as f:
            lines= f.readlines()
        with open('exAA.yaml', 'w') as f:
            for line in lines:
                if 'title:' in line:
                    f.write('title: exAA_updated\n')
                else:
                    f.write(line)
        self.manager.update_example(1,'exAA')
        os.chdir('..')
        shutil.rmtree('exa_updated')
        with open(os.path.join(self.manager.path, 'info.yaml'), 'r') as f:
            info = yaml.safe_load(f)
        self.assertIn('examples', info)
        self.assertEqual(len(info['examples']), 2)
        self.assertEqual(info['examples'][0]['name'], 'exAA')
        self.assertEqual(info['examples'][0]['index'], 1)
        os.chdir(cwd)
        