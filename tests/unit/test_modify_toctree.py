import unittest
from pestifer.sphinxext.modify_toctree import modify_toctree

class TestModifyToctree(unittest.TestCase):
    def setUp(self):
        self.sample_lines = [
            ".. toctree::\n",
            "  :maxdepth: 2\n",
            "\n",
            "  howlers/example1\n",
            "  howlers/example2\n",
            "  howlers/example3\n",
            "\n",
        ]
        with open('sample_rst_file.rst', 'w') as f:
            f.writelines(self.sample_lines)

    def test_modify_toctree_add(self):
        filepath = 'sample_rst_file.rst'
        action = "add"
        new_entry = "example4"
        modify_toctree(filepath, action, new_entry=new_entry)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        self.assertIn("  howlers/example4\n", lines)

    def test_modify_toctree_delete(self):
        filepath = 'sample_rst_file.rst'
        action = "delete"
        target = "howlers/example1"
        modify_toctree(filepath, action, target=target)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        self.assertNotIn("  example1\n", lines)

    def test_modify_toctree_insert(self):
        filepath = 'sample_rst_file.rst'
        action = "insert"
        target = "howlers/example1"
        new_entry = "example2.5"
        modify_toctree(filepath, action, target=target, new_entry=new_entry)
        with open(filepath, 'r') as f:
            lines = f.readlines()
        self.assertIn("  howlers/example2.5\n", lines)
        idx = lines.index("  howlers/example2.5\n")
        self.assertEqual(lines[idx - 1].strip(), "howlers/example1")