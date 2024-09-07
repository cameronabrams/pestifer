# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import unittest
from pestifer.command import Command

class TestCommand(unittest.TestCase):
    def test_echo(self):
        expected_stdout='123\n'
        c=Command('echo 123')
        c.run()
        self.assertEqual(expected_stdout,c.stdout)
