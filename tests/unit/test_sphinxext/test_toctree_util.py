import unittest
from pestifer.sphinxext.toctree_util import modify_toctree, get_name_from_toctree, get_num_entries_in_toctree

import os
import shutil
class TestModifyToctree(unittest.TestCase):
    def setUp(self):
        shutil.copy('../fixtures/toctree_utils_inputs/sample_rst_file.rst', '.')
        shutil.copy('../fixtures/toctree_utils_inputs/empty_rst_file.rst', '.')
        self.rst_folder_path = 'examples'

    def tearDown(self):
        if os.path.isfile('sample_rst_file.rst'):
            os.remove('sample_rst_file.rst')
        if os.path.isfile('empty_rst_file.rst'):
            os.remove('empty_rst_file.rst')

    def test_get_num_entries_in_toctree(self):
        filepath = 'sample_rst_file.rst'
        num_entries = get_num_entries_in_toctree(filepath)
        self.assertEqual(num_entries, 3)

        filepath = 'empty_rst_file.rst'
        num_entries = get_num_entries_in_toctree(filepath)
        self.assertEqual(num_entries, 0)

    def test_new_entry_new_toctree(self):
        filepath = 'empty_rst_file.rst'
        action = "append"
        new_entry = "ex01/exampleD"
        modify_toctree(filepath, action, new_entry=new_entry, common_prefix=self.rst_folder_path)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        self.assertIn(f"  {self.rst_folder_path}/ex01/exampleD\n", lines)
        self.assertEqual(get_num_entries_in_toctree(filepath), 1)

    def test_modify_toctree_add_entry(self):
        filepath = 'sample_rst_file.rst'
        action = "append"
        new_entry = "ex04/example4"
        modify_toctree(filepath, action, new_entry=new_entry)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        self.assertIn("  howlers/ex04/example4\n", lines[-1])

    def test_modify_toctree_delete_entry(self):
        filepath = 'sample_rst_file.rst'
        action = "delete"
        target = "ex02/exampleB"
        modify_toctree(filepath, action, target=target)
        with open(filepath, 'r') as f:
            all_content = f.read()
        self.assertNotIn("  howlers/ex02/exampleB\n", all_content)

    def test_modify_toctree_update_entry(self):
        filepath = 'sample_rst_file.rst'
        action = "update"
        target = "ex01/exampleA"
        new_entry = "ex01/exampleA_updated"
        modify_toctree(filepath, action, target=target, new_entry=new_entry)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        self.assertIn("  howlers/ex01/exampleA_updated\n", lines)
        self.assertNotIn("  howlers/ex01/exampleA\n", lines)

    def test_get_name_from_toctree(self):
        filepath = 'sample_rst_file.rst'
        match_str = "ex02"
        name = get_name_from_toctree(filepath, match_str)
        self.assertEqual(name, "exampleB")