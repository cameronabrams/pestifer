import unittest
import shutil
import os
from pathlib import Path
from pestifer.core.examplemanager import ExampleManager

class TestExampleManager(unittest.TestCase):

    def setUp(self):
        assert os.getcwd().endswith('tests/unit/test_core/test_examplemanager')
        assert os.path.isdir('../fixtures')
        assert os.path.isdir('../fixtures/test_examplemanager_inputs')
        assert os.path.isdir('../fixtures/test_examplemanager_inputs/project')
        assert os.path.isdir('../fixtures/test_examplemanager_inputs/userspace')
        if os.path.isdir('./project'):
            shutil.rmtree('./project', ignore_errors=True)
        if os.path.isdir('./userspace'):
            shutil.rmtree('./userspace', ignore_errors=True)
        shutil.copytree('../fixtures/test_examplemanager_inputs/project', './project/', dirs_exist_ok=True)
        shutil.copytree('../fixtures/test_examplemanager_inputs/userspace', './userspace/', dirs_exist_ok=True)
        self.manager = ExampleManager(
            examples_path='project/package/resources/examples',
            docs_source_path='project/docs/source'
        )

    def tearDown(self):
        shutil.rmtree('project')
        shutil.rmtree('userspace')

    def test_example_manager_init(self):
        self.assertIsNotNone(self.manager)
        self.assertEqual(len(self.manager.examples), 0)
        self.assertEqual(self.manager.path, Path(os.path.join(os.getcwd(), 'project', 'package', 'resources', 'examples')))
        self.assertTrue(os.path.isdir(self.manager.path))

    def _build_example_set(self):
        os.chdir('userspace')
        self.example1 = self.manager.append_example(1, 'exA.yaml', auxiliary_inputs=['exA_companion1.pdb', 'exA_companion1.psf'])
        self.example2 = self.manager.append_example(2, 'exB.yaml')
        os.chdir('..')

    def test_example_manager_build_example_set(self):
        self._build_example_set()
        self.assertEqual(len(self.manager.examples), 2)
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[0].scriptpath)))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[0].inputspath, 'exA_companion1.pdb')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[0].inputspath, 'exA_companion1.psf')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[1].inputspath, 'exB.yaml')))
        self.assertTrue(os.path.isfile(self.manager.sphinx_example_manager.examples_rst))
        self.assertTrue(os.path.isdir(self.manager.sphinx_example_manager.examples_folder_path))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.sphinx_example_manager.examples_folder_path, 'ex01', 'exA.rst')))
        self.assertTrue(os.path.isfile(os.path.join(self.manager.sphinx_example_manager.examples_folder_path, 'ex02', 'exB.rst')))
    
    def test_example_manager_delete_example(self):
        self._build_example_set()
        self.manager.delete_example(2)
        self.assertEqual(len(self.manager.examples), 1)
        self.assertTrue(os.path.isfile(os.path.join(self.manager.path, self.manager.examples[0].scriptpath)))
        self.assertFalse(os.path.isfile(os.path.join(self.manager.path, self.example2.scriptpath)))
        self.assertFalse(os.path.isfile(os.path.join(self.manager.sphinx_example_manager.examples_folder_path, 'ex02', 'exB.rst')))

    def test_example_manager_update_example_inplace(self):
        self._build_example_set()
        ex_rst = os.path.join(self.manager.sphinx_example_manager.examples_folder_path, 'ex01', 'exA.rst')
        self.assertTrue(os.path.isfile(ex_rst))
        with open(ex_rst, 'r') as f:
            rst_current = f.read()
        self.assertNotIn('PDB ID 9ggs', rst_current)
        os.chdir('userspace')
        os.mkdir('exa_updated')
        os.chdir('exa_updated')
        Path('exA_companion1.pdb').touch()
        Path('exA_companion1.psf').touch()
        self.manager.update_example(1, title='Updated title', db_id='9ggs', author_name='Mel Brooks')
        os.chdir('..')
        shutil.rmtree('exa_updated')
        os.chdir('..')
        self.assertTrue(os.path.exists('project'))
        example = self.manager.examples.get_example_by_example_id(1)
        self.assertEqual(example.title, 'Updated title')
        self.assertEqual(example.author_name, 'Mel Brooks')
        self.assertEqual(example.shortname, 'exA')
        ex_rst = os.path.join(self.manager.sphinx_example_manager.examples_folder_path, 'ex01', 'exA.rst')
        self.assertTrue(os.path.isfile(ex_rst))
        with open(ex_rst, 'r') as f:
            rst_current = f.read()
        self.assertIn('PDB ID 9ggs', rst_current)

    def test_example_manager_update_example_rename(self):
        cwd = os.getcwd()
        self._build_example_set()
        os.chdir('userspace')
        os.mkdir('exa_updated')
        os.chdir('exa_updated')
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
        os.chdir(cwd)
        example = self.manager.examples.get_example_by_example_id(1)
        self.assertEqual(example.title, 'exAA_updated')
        self.assertEqual(example.shortname, 'exAA')
        ex_rst = os.path.join(self.manager.sphinx_example_manager.examples_folder_path, 'ex01', 'exAA.rst')
        self.assertTrue(os.path.isfile(ex_rst))
        with open(ex_rst, 'r') as f:
            rst_current = f.read()
        self.assertIn('examples/01/inputs/exAA.yaml', rst_current)
        