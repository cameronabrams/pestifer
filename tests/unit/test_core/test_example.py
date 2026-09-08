import unittest
from pathlib import Path
from pestifer.core.example import Example, ExampleList
from pestifer.core.errors import PestiferError

class TestExample(unittest.TestCase):
    
    def test_example_creation(self):
        example = Example(example_id=1, shortname="example1", db_id='abc123', title="This is an example")
        self.assertEqual(example.example_id, 1)
        self.assertEqual(example.shortname, "example1")
        self.assertEqual(example.db_id, 'abc123')
        self.assertEqual(example.title, "This is an example")
        self.assertEqual(example.inputspath, Path('01/inputs'))
        self.assertEqual(example.outputspath, Path('01/outputs'))
        self.assertEqual(example.scriptpath, Path('01/inputs/example1.yaml'))

    def test_example_report_line(self):
        example = Example(example_id=1, shortname="example1", db_id='abc1', title="This is an example")
        report_line = example.report_line()
        self.assertEqual(report_line, "      1    abc1  example1                          This is an example")

    def test_example_list_creation(self):
        example1 = Example(example_id=1, shortname="example1", db_id='abc123', title="This is an example")
        example2 = Example(example_id=2, shortname="example2", db_id='def456', title="This is another example")
        example3 = Example(example_id=3, shortname="example3", db_id='ghi789', title="This is yet another example")
        example_list = ExampleList([example1, example2])
        self.assertEqual(len(example_list), 2)
        self.assertEqual(example_list[0].shortname, "example1")
        self.assertEqual(example_list[1].shortname, "example2")
        self.assertIn(example1, example_list)
        self.assertIn(example2, example_list)
        example_list.append(example3)
        self.assertEqual(len(example_list), 3)
        self.assertEqual(example_list[2].shortname, "example3")
        self.assertIn(example3, example_list)
    
    def test_example_list_insert(self):
        example1 = Example(example_id=1, shortname="example1", db_id='abc123', title="This is an example")
        example2 = Example(example_id=2, shortname="example2", db_id='def456', title="This is another example")
        example3 = Example(example_id=3, shortname="example3", db_id='ghi789', title="This is yet another example")
        example_list = ExampleList([example1, example2])
        example_list.insert(1, example3)
        self.assertEqual(len(example_list), 3)
        self.assertEqual(example_list[1].shortname, "example3")
        self.assertEqual(example_list[0].shortname, "example1")
        self.assertEqual(example_list[2].shortname, "example2")

    def test_example_list_remove(self):
        example1 = Example(example_id=1, shortname="example1", db_id='abc123', title="This is an example")
        example2 = Example(example_id=2, shortname="example2", db_id='def456', title="This is another example")
        example3 = Example(example_id=3, shortname="example3", db_id='ghi789', title="This is yet another example")
        example_list = ExampleList([example1, example2, example3])
        self.assertEqual(len(example_list), 3)
        example_list.remove(example2)
        self.assertEqual(len(example_list), 2)
        self.assertNotIn(example2, example_list)

    def test_example_list_read_yaml(self):
        example1 = Example(example_id=1, shortname="example1", db_id='abc123', title="This is an example")
        example2 = Example(example_id=2, shortname="example2", db_id='def456', title="This is another example")
        example_list = ExampleList([example1, example2])
        yaml_content = example_list.to_yaml()
        loaded_list = ExampleList.read_yaml(yaml_content)
        self.assertEqual(len(loaded_list), 2)
        self.assertEqual(loaded_list[0].shortname, "example1")
        self.assertEqual(loaded_list[1].shortname, "example2")

class TestExampleResolve(unittest.TestCase):
    """
    ``resolve`` lets a caller name an example by id or by shortname.  The two orderings the
    examples are presented in -- by biological subject in the repo, by pestifer capability in
    the paper -- cannot both be made sequential by one integer, so the shortname is the handle
    that does not depend on either.
    """

    def setUp(self):
        self.examples = ExampleList([
            Example(example_id=23, shortname='subtilisin-dmso', title='dmso'),
            Example(example_id=24, shortname='subtilisin-acetone', title='acetone'),
            Example(example_id=26, shortname='subtilisin-acetonitrile', title='acetonitrile'),
        ])

    def test_resolves_int_id(self):
        self.assertEqual(self.examples.resolve(24).shortname, 'subtilisin-acetone')

    def test_resolves_digit_string_as_id(self):
        # argparse hands the token over as a string; '24' must stay an id, not become a name
        self.assertEqual(self.examples.resolve('24').shortname, 'subtilisin-acetone')

    def test_resolves_shortname(self):
        self.assertEqual(self.examples.resolve('subtilisin-acetone').example_id, 24)

    def test_unknown_id_reports_the_range(self):
        with self.assertRaises(PestiferError) as cm:
            self.examples.resolve(99)
        self.assertIn('ids run 1-26', str(cm.exception))

    def test_unknown_name_suggests_near_misses(self):
        with self.assertRaises(PestiferError) as cm:
            self.examples.resolve('subtilisin-aceton')
        msg = str(cm.exception)
        self.assertIn('did you mean', msg)
        self.assertIn('subtilisin-acetone', msg)

    def test_unrecognizable_name_points_at_the_listing(self):
        with self.assertRaises(PestiferError) as cm:
            self.examples.resolve('zzzzzzzz')
        self.assertIn('show-resources examples', str(cm.exception))

    def test_raises_pestifer_error_not_keyerror(self):
        # the CLI turns PestiferError into a one-line message; a KeyError would traceback
        with self.assertRaises(PestiferError):
            self.examples.resolve('nope')
