# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Exit status derived from a subcommand's return value.

Most subcommands return ``True`` to mean "done", and ``bool`` is a subclass of ``int``, so a
naive int check turns every one of them into a failing exit code -- which is exactly what
happened to ``pestifer show-resources``.  Only a genuine int, or a genuine int ``exit_code``
attribute, is a status.
"""
import unittest

from pestifer.cli.pestifer import subcommand_exit_code


class _Result:
    def __init__(self, exit_code):
        self.exit_code = exit_code


class TestSubcommandExitCode(unittest.TestCase):

    def test_success_sentinels_are_zero(self):
        for value in (None, True, False, 0):
            with self.subTest(value=value):
                self.assertEqual(subcommand_exit_code(value), 0)

    def test_true_is_not_an_exit_code(self):
        # bool is a subclass of int; returning True must not mean "exit 1"
        self.assertEqual(subcommand_exit_code(True), 0)

    def test_int_return_is_the_status(self):
        self.assertEqual(subcommand_exit_code(1), 1)
        self.assertEqual(subcommand_exit_code(127), 127)

    def test_exit_code_attribute_is_honored(self):
        self.assertEqual(subcommand_exit_code(_Result(0)), 0)
        self.assertEqual(subcommand_exit_code(_Result(3)), 3)

    def test_non_integer_results_are_success(self):
        self.assertEqual(subcommand_exit_code('text'), 0)
        self.assertEqual(subcommand_exit_code(_Result(True)), 0)
        self.assertEqual(subcommand_exit_code(object()), 0)


if __name__ == '__main__':
    unittest.main()
