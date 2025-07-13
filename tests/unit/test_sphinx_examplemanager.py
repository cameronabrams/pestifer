import unittest
import os
import shutil
from pestifer.sphinxext.sphinx_examplemanager import SphinxExampleManager
from pestifer.sphinxext.toctree_util import modify_toctree, get_num_entries_in_toctree
from pestifer.core.example import Example

class TestSphinxExampleManager(unittest.TestCase):
    def setUp(self):
        if os.path.isdir('docs'):
            shutil.rmtree('docs', ignore_errors=True)
        if os.path.isdir('userspace'):
            shutil.rmtree('userspace', ignore_errors=True)
        shutil.copytree('../fixtures/sphinx_example_inputs', '.', dirs_exist_ok=True)
        self.manager = SphinxExampleManager(docs_source_path='docs/source',examples_folder_name='examples', examples_rst_name='examples.rst')

    def _build_example_set(self):
        self.examples=[]
        for l in ['A','B','C']:
            self.examples.append(Example.from_yaml(f'userspace/ex{l}.yaml'))
            self.manager.insert_example(len(self.examples), self.examples[-1])
        self.updated_existing_example = Example.from_yaml('userspace/exA_updated.yaml')
        self.example_to_insert=Example.from_yaml('userspace/exD.yaml')

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
        rst_files= os.listdir(self.manager.examples_folder_path)
        self.assertTrue(len(rst_files) == 0)
        self._build_example_set()
        self.assertTrue(get_num_entries_in_toctree(self.manager.examples_rst) == 3)
        rst_files= os.listdir(self.manager.examples_folder_path)
        self.assertIn('exA.rst', rst_files)
        self.assertIn('exB.rst', rst_files)
        self.assertIn('exC.rst', rst_files)
        self.assertTrue(len(rst_files) == 3)

    def test_sphinx_example_manager_insert_example(self):
        self._build_example_set()
        self.example_to_insert.index = 2
        self.manager.insert_example(2,self.example_to_insert)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            toctree_lines = f.readlines()
        exAidx=toctree_lines.index('  examples/exA\n')
        exINSidx=toctree_lines.index('  examples/exD\n')
        exBidx=toctree_lines.index('  examples/exB\n')
        self.assertTrue(exAidx < exINSidx < exBidx)
        rst_files= os.listdir(self.manager.examples_folder_path)
        self.assertIn('exD.rst', rst_files)

    def test_sphinx_example_manager_append_example(self):
        self._build_example_set()
        self.manager.append_example(self.example_to_insert)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            toctree_lines = f.readlines()
        exCidx=toctree_lines.index('  examples/exC\n')
        exINSidx=toctree_lines.index('  examples/exD\n')
        self.assertTrue(exCidx < exINSidx)
        rst_files= os.listdir(self.manager.examples_folder_path)
        self.assertIn('exD.rst', rst_files)

    def test_sphinx_example_manager_delete_example(self):
        self._build_example_set()
        example_to_delete = self.examples[1]
        self.assertTrue(example_to_delete.name == 'exB')
        self.manager.delete_example(example_to_delete)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            toctree_lines = f.readlines()
        self.assertNotIn('  examples/exB\n', toctree_lines)
        rst_files= os.listdir(self.manager.examples_folder_path)
        self.assertNotIn('exB.rst', rst_files)
         
    def test_sphinx_example_manager_title_index_shift(self):
        self._build_example_set()
        self.manager.title_index_shift(start_index=2, shift=1)
        example_rst_path = os.path.join(self.manager.examples_folder_path, f'exB.rst')
        self.assertTrue(os.path.isfile(example_rst_path))
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 3:", content)
        self.manager.title_index_shift(start_index=2, shift=-1)
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 2:", content)

    def test_sphinx_example_manager_append_existing_example(self):
        self._build_example_set()
        with self.assertRaises(FileExistsError):
            self.manager.append_example(self.examples[0])  # Attempt to add an existing example
    
    def test_sphinx_example_manager_update_example(self):
        self._build_example_set()
        self.manager.update_example(1, self.updated_existing_example)
        example_rst_path = os.path.join(self.manager.examples_folder_path, f'{self.updated_existing_example.name}.rst')
        self.assertTrue(os.path.isfile(example_rst_path))
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("PDB ID 1ccc", content)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            content = f.read()
        self.assertIn(f'examples/{self.updated_existing_example.name}', content)
        content=content.replace(f'examples/{self.updated_existing_example.name}', 'PLACEHOLDER')
        self.assertNotIn(f'examples/{self.updated_existing_example.name}', content)