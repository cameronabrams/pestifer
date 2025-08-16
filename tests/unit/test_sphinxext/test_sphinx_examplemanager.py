import os
import shutil
import unittest

from pestifer.sphinxext.sphinx_examplemanager import SphinxExampleManager
from pestifer.sphinxext.toctree_util import get_num_entries_in_toctree
from pestifer.core.example import Example

class TestSphinxExampleManager(unittest.TestCase):

    def setUp(self):
        if os.path.isdir('project'):
            shutil.rmtree('project', ignore_errors=True)
        if os.path.isdir('userspace'):
            shutil.rmtree('userspace', ignore_errors=True)
        shutil.copytree('../fixtures/sphinx_example_inputs', '.', dirs_exist_ok=True)
        self.manager = SphinxExampleManager(docs_source_path='project/docs/source', examples_folder_name='examples', examples_rst_name='examples.rst')

    def _build_example_set(self):
        self.examples=[]
        for i,l in enumerate(['A','B','C']):
            title, db_id = Example._get_implied_metadata(f'userspace/ex{l}.yaml')
            example = Example(example_id=i+1, title=title, db_id=db_id, shortname=f'ex{l}')
            self.examples.append(example)
            self.manager.append_example(example)
        title, db_id = Example._get_implied_metadata('userspace/exA_updated.yaml')
        self.updated_existing_example = Example(example_id=1, title=title, db_id=db_id, shortname='exA_updated')
        title, db_id = Example._get_implied_metadata('userspace/exD.yaml')
        self.example_to_add = Example(example_id=4, title=title, db_id=db_id, shortname='exD')

    def tearDown(self):
        shutil.rmtree('docs', ignore_errors=True)
        shutil.rmtree('userspace', ignore_errors=True)
        return super().tearDown()

    def test_sphinx_example_manager_init(self):
        self.assertIsNotNone(self.manager)
        self.assertTrue(os.path.isdir(self.manager.docs_source_path))
        self.assertTrue(os.path.isfile(self.manager.examples_rst))
        self.assertTrue(os.path.isdir(self.manager.examples_folder_path))
        self.assertTrue(get_num_entries_in_toctree(self.manager.examples_rst) == 0)
        ex_folders = os.listdir(self.manager.examples_folder_path)
        self.assertTrue(len(ex_folders) == 0)
        self._build_example_set()
        self.assertTrue(get_num_entries_in_toctree(self.manager.examples_rst) == 3)
        ex_folders = os.listdir(self.manager.examples_folder_path)
        self.assertIn('ex01', ex_folders)
        self.assertIn('ex02', ex_folders)
        self.assertIn('ex03', ex_folders)
        self.assertTrue(len(ex_folders) == 3)

    def test_sphinx_example_manager_append_example(self):
        self._build_example_set()
        self.manager.append_example(self.example_to_add)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            toctree_lines = f.readlines()
        exCidx=toctree_lines.index('  examples/ex03/exC\n')
        exINSidx=toctree_lines.index('  examples/ex04/exD\n')
        self.assertTrue(exCidx < exINSidx)
        example_folders = os.listdir(self.manager.examples_folder_path)
        self.assertIn('ex04', example_folders)
        self.assertTrue(len(example_folders) == 4)
        new_example_folder = os.path.join(self.manager.examples_folder_path, 'ex04')
        self.assertTrue(os.path.exists(os.path.join(new_example_folder, 'exD.rst')))

    def test_sphinx_example_manager_delete_example(self):
        self._build_example_set()
        example_to_delete = self.examples[1]
        self.assertTrue(example_to_delete.shortname == 'exB')
        self.manager.delete_example(example_to_delete)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            toctree_lines = f.readlines()
        self.assertNotIn(f'  examples/{example_to_delete.rootfolderpath}/exB\n', toctree_lines)
        example_folders = os.listdir(self.manager.examples_folder_path)
        self.assertNotIn(example_to_delete.rootfolderpath, example_folders)
        self.assertEqual(len(example_folders), 2)

    def test_sphinx_example_manager_append_existing_example(self):
        self._build_example_set()
        with self.assertRaises(AssertionError):
            self.manager.append_example(self.examples[0])  # Attempt to add an existing example
    
    def test_sphinx_example_manager_update_example(self):
        self._build_example_set()
        self.manager.update_example(1, self.updated_existing_example)
        example_rst_path = os.path.join(self.manager.examples_folder_path, self.updated_existing_example.rootfolderpath, f'{self.updated_existing_example.shortname}.rst')
        self.assertTrue(os.path.isfile(example_rst_path))
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("PDB ID 1ccc", content)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            content = f.read()
        self.assertIn(f'examples/{self.updated_existing_example.rootfolderpath}/{self.updated_existing_example.shortname}', content)
        content=content.replace(f'examples/{self.updated_existing_example.rootfolderpath}/{self.updated_existing_example.shortname}', 'PLACEHOLDER')
        self.assertNotIn(f'examples/{self.updated_existing_example.rootfolderpath}/{self.updated_existing_example.shortname}', content)