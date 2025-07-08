import unittest
from pestifer.sphinxext.toctree_util import modify_toctree,get_name_from_toctree

import os
import shutil
class TestModifyToctree(unittest.TestCase):
    def setUp(self):
        shutil.copy('../fixtures/toctree_utils_inputs/sample_rst_file.rst', '.')

    def tearDown(self):
        if os.path.isfile('sample_rst_file.rst'):
            os.remove('sample_rst_file.rst')

    def test_modify_toctree_add(self):
        filepath = 'sample_rst_file.rst'
        action = "add"
        new_entry = "example4"
        modify_toctree(filepath, action, new_entry=new_entry)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        self.assertIn("  howlers/example4\n", lines[-1])

    def test_modify_toctree_delete(self):
        filepath = 'sample_rst_file.rst'
        action = "delete"
        target = "howlers/exampleB"
        modify_toctree(filepath, action, target=target)
        with open(filepath, 'r') as f:
            all_content = f.read()
        self.assertNotIn("  exampleB\n", all_content)

    def test_modify_toctree_insert(self):
        filepath = 'sample_rst_file.rst'
        action = "insert"
        new_index=2 # 1-based index for insertion
        new_entry = "exampleD"
        modify_toctree(filepath, action, new_entry=new_entry, new_index=new_index)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        self.assertIn("  howlers/exampleD\n", lines)
        idx = lines.index("  howlers/exampleD\n")
        self.assertEqual(lines[idx - 1].strip(), "howlers/exampleA")
        self.assertEqual(lines[idx + 1].strip(), "howlers/exampleB")

    def test_modify_toctree_update(self):
        filepath = 'sample_rst_file.rst'
        action = "update"
        target = "howlers/exampleA"
        new_entry = "exampleA_updated"
        modify_toctree(filepath, action, target=target, new_entry=new_entry)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        self.assertIn("  howlers/exampleA_updated\n", lines)
        self.assertNotIn("  howlers/exampleA\n", lines)

    def test_get_name_from_toctree(self):
        filepath = 'sample_rst_file.rst'
        index = 2
        name = get_name_from_toctree(filepath, index)
        self.assertEqual(name, "exampleB")