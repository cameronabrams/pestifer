import unittest
import os
import yaml
from pestifer.core.example import Example, ExampleList

class TestExample(unittest.TestCase):
    def test_example_creation(self):
        example = Example(name="example1", pdbID='abc123', description="This is an example", index=1)
        self.assertEqual(example.name, "example1")
        self.assertEqual(example.pdbID, 'abc123')
        self.assertEqual(example.description, "This is an example")
        self.assertEqual(example.index, 1)
        with open(example.name+'.yaml', 'w') as f:
            f.write(example.to_yaml())
        self.assertTrue(os.path.exists(example.name+'.yaml'))
        with open(example.name+'.yaml', 'r') as f:
            content = yaml.safe_load(f)
        self.assertEqual(content['name'], "example1")
        self.assertEqual(content['pdbID'], 'abc123')
        self.assertEqual(content['description'], "This is an example")
        self.assertEqual(content['index'], 1)
        os.remove(example.name+'.yaml')

    def test_example_report_line(self):
        example = Example(name="example1", pdbID='abc1', description="This is an example", index=1)
        report_line = example.report_line()
        self.assertEqual(report_line, "      1    abc1  example1                          This is an example")

    def test_example_list_creation(self):
        example1 = Example(name="example1", pdbID='abc123', description="This is an example", index=1)
        example2 = Example(name="example2", pdbID='def456', description="This is another example", index=2)
        example3 = Example(name="example3", pdbID='ghi789', description="This is yet another example", index=3)
        example_list = ExampleList([example1, example2])
        self.assertEqual(len(example_list), 2)
        self.assertEqual(example_list[0].name, "example1")
        self.assertEqual(example_list[1].name, "example2")
        self.assertIn(example1, example_list)
        self.assertIn(example2, example_list)
        example_list.append(example3)
        self.assertEqual(len(example_list), 3)
        self.assertEqual(example_list[2].name, "example3")
        self.assertIn(example3, example_list)
    
    def test_example_list_insert(self):
        example1 = Example(name="example1", pdbID='abc123', description="This is an example", index=1)
        example2 = Example(name="example2", pdbID='def456', description="This is another example", index=2)
        example3 = Example(name="example3", pdbID='ghi789', description="This is yet another example", index=3)
        example_list = ExampleList([example1, example2])
        example_list.insert(1, example3)
        self.assertEqual(len(example_list), 3)
        self.assertEqual(example_list[1].name, "example3")
        self.assertEqual(example_list[0].name, "example1")
        self.assertEqual(example_list[2].name, "example2")

    def test_example_list_remove(self):
        example1 = Example(name="example1", pdbID='abc123', description="This is an example", index=1)
        example2 = Example(name="example2", pdbID='def456', description="This is another example", index=2)
        example3 = Example(name="example3", pdbID='ghi789', description="This is yet another example", index=3)
        example_list = ExampleList([example1, example2, example3])
        self.assertEqual(len(example_list), 3)
        example_list.remove(example2)
        self.assertEqual(len(example_list), 2)
        self.assertNotIn(example2, example_list)

    def test_example_list_read_yaml(self):
        example1 = Example(name="example1", pdbID='abc123', description="This is an example", index=1)
        example2 = Example(name="example2", pdbID='def456', description="This is another example", index=2)
        example_list = ExampleList([example1, example2])
        yaml_content = example_list.to_yaml()
        loaded_list = ExampleList.read_yaml(yaml_content)
        self.assertEqual(len(loaded_list), 2)
        self.assertEqual(loaded_list[0].name, "example1")
        self.assertEqual(loaded_list[1].name, "example2")