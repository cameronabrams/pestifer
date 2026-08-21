# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the ``catdcd`` version gate.

Deliberately *not* marked ``needs_tools``: the whole point of the check is to catch a toolchain
that is present but wrong, so the logic must be testable where no toolchain exists at all --
which is exactly where CI runs.  ``Config`` is built with ``__new__`` and given only the one
attribute the check reads.
"""

import unittest
from unittest import mock

from pestifer.core.config import Config, MIN_CATDCD_VERSION, _version_tuple
from pestifer.core.errors import PestiferError


class TestVersionTuple(unittest.TestCase):

    def test_parses_dotted_versions(self):
        self.assertEqual(_version_tuple('5.2'), (5, 2))
        self.assertEqual(_version_tuple('  5.10 '), (5, 10))
        self.assertEqual(_version_tuple('5.2beta'), (5, 2))

    def test_bare_major_is_minor_zero(self):
        self.assertEqual(_version_tuple('5'), (5, 0))

    def test_unparseable_is_none(self):
        for text in ('unknown', '', None, 'CatDCD'):
            self.assertIsNone(_version_tuple(text))


class TestCatdcdVersionGate(unittest.TestCase):

    def _config(self, cmd='catdcd'):
        c = Config.__new__(Config)
        c.shell_commands = {'catdcd': cmd}
        return c

    def _check(self, reported, resolves=True):
        with mock.patch('pestifer.util.provenance.catdcd_version', return_value=reported), \
             mock.patch('pestifer.core.config.shutil.which',
                        return_value='/usr/local/bin/catdcd' if resolves else None):
            return self._config()._verify_catdcd_version()

    def test_current_version_passes(self):
        self.assertIsNone(self._check('5.2'))

    def test_newer_version_passes(self):
        # guards against comparing version *strings*: '5.10' < '5.2' lexically, but 5.10 is newer
        self.assertIsNone(self._check('5.10'))

    def test_old_version_is_refused(self):
        with self.assertRaises(PestiferError) as cm:
            self._check('5.1')
        msg = str(cm.exception)
        self.assertIn('5.1', msg)
        self.assertIn('insertion code', msg)   # say *why*, not just 'too old'

    def test_much_older_version_is_refused(self):
        with self.assertRaises(PestiferError):
            self._check('4.9')

    def test_unknown_version_warns_but_continues(self):
        # an advisory probe must never block a build through its own malfunction
        with self.assertLogs('pestifer.core.config', level='WARNING') as cm:
            self.assertIsNone(self._check('unknown'))
        self.assertIn('insertion code', ''.join(cm.output))

    def test_unresolvable_command_is_not_this_checks_business(self):
        # presence is enforced by the required-command loop, which raises its own error
        self.assertIsNone(self._check('5.1', resolves=False))

    def test_minimum_is_the_documented_5_2(self):
        # docs/source/installation.rst states 5.2 as the requirement; keep them in step
        self.assertEqual(MIN_CATDCD_VERSION, (5, 2))


if __name__ == '__main__':
    unittest.main()
