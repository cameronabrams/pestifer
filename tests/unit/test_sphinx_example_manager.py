import unittest
import os
import shutil
import tarfile
from pestifer.sphinxext.sphinx_examplemanager import SphinxExampleManager
from pestifer.core.example import Example

class TestSphinxExampleManager(unittest.TestCase):
    def setUp(self):
        if os.path.isdir('docs'):
            shutil.rmtree('docs', ignore_errors=True)
        if os.path.isfile('example1.yaml'):
            os.remove('example1.yaml')
        if os.path.isfile('example2-5.yaml'):
            os.remove('example2-5.yaml')
        shutil.copytree('../fixtures/sphinx_example_inputs', '.', dirs_exist_ok=True)
        self.manager = SphinxExampleManager(docs_source_path='docs/source')
        self.example = Example(name="example2-5", pdbID='abc123', description="This is an example")
        self.existing_example=Example(name="example1", pdbID='def456', description="Example Input (orig 1)", index=1)
        self.updated_existing_example=Example(name="example1_updated", pdbID='def457', description="Example Input (orig 1) updated", index=1)

    def tearDown(self):
        shutil.rmtree('docs', ignore_errors=True)
        os.remove('example1.yaml')
        os.remove('example2-5.yaml')
        return super().tearDown()

    def test_sphinx_example_manager_init(self):
        self.assertIsNotNone(self.manager)
        self.assertTrue(os.path.isdir(self.manager.docs_source_path))
        self.assertTrue(os.path.isfile(self.manager.examples_rst))
        self.assertTrue(os.path.isdir(self.manager.examples_source_path))
        self.assertTrue(os.path.exists('example2-5.yaml'))

    def test_sphinx_example_manager_add_example(self):
        self.example.index=4 # this will be set by the manager
        self.manager.add_example(self.example)
        example_rst_path = os.path.join(self.manager.examples_source_path, f'{self.example.name}.rst')
        self.assertTrue(os.path.isfile(example_rst_path))
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 4: This is an example", content)
        self.assertIn("PDB ID abc123", content)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            toctree_lines = f.readlines()
        self.assertIn(f'examples/example2-5', toctree_lines[-1])

    def test_sphinx_example_manager_title_index_shift(self):
        self.manager.title_index_shift(start_index=2, shift=1)
        example_rst_path = os.path.join(self.manager.examples_source_path, f'example2.rst')
        self.assertTrue(os.path.isfile(example_rst_path))
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 3: Example Input (orig 2)", content)
        self.manager.title_index_shift(start_index=2, shift=-1)
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 2: Example Input (orig 2)", content)

    def test_sphinx_example_manager_insert_delete_example(self):
        self.example.index = 2
        self.manager.insert_example(2,self.example)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            toctree_lines = f.readlines()
        ex1idx=toctree_lines.index('  examples/example1\n')
        exINSidx=toctree_lines.index('  examples/example2-5\n')
        ex2idx=toctree_lines.index('  examples/example2\n')
        self.assertTrue(ex1idx < exINSidx < ex2idx)
        example_rst_path = os.path.join(self.manager.examples_source_path, f'{self.example.name}.rst')
        self.assertTrue(os.path.isfile(example_rst_path))
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 2: This is an example", content)
        self.assertIn("PDB ID abc123", content)
        existing_example_rst_path = os.path.join(self.manager.examples_source_path, f'example2.rst')
        self.assertTrue(os.path.isfile(existing_example_rst_path))
        with open(existing_example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 3: Example Input (orig 2)", content)
        self.manager.delete_example(self.example)
        self.assertFalse(os.path.isfile(example_rst_path))
        with open(existing_example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 2: Example Input (orig 2)", content)

    def test_sphinx_example_manager_add_existing_example(self):
        with self.assertRaises(FileExistsError):
            self.manager.add_example(self.existing_example)
        # Check that the example was not added again
        example_rst_path = os.path.join(self.manager.examples_source_path, f'{self.existing_example.name}.rst')
        self.assertTrue(os.path.isfile(example_rst_path))
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 1:", content)
    
    def test_sphinx_example_manager_update_example(self):
        self.manager.update_example(1, self.updated_existing_example)
        example_rst_path = os.path.join(self.manager.examples_source_path, f'{self.updated_existing_example.name}.rst')
        self.assertTrue(os.path.isfile(example_rst_path))
        with open(example_rst_path, 'r') as f:
            content = f.read()
        self.assertIn("Example 1: Example Input (orig 1)", content)
        self.assertIn("PDB ID def457", content)
        toctree_path = self.manager.examples_rst
        with open(toctree_path, 'r') as f:
            content = f.read()
        self.assertIn(f'examples/{self.updated_existing_example.name}', content)
        content=content.replace(f'examples/{self.updated_existing_example.name}', 'PLACEHOLDER')
        self.assertNotIn(f'examples/{self.existing_example.name}', content)