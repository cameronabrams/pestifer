# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the configurable Open Babel executable.

Every other external command pestifer invokes -- ``vmd``, ``namd3``, ``charmrun``, ``catdcd`` --
can be relocated through a ``paths`` key, for the user whose copy lives under another name or
outside ``PATH`` (a module-loaded HPC build, a conda environment).  ``obabel`` could not, because
``mol2_writer`` hardcoded the string.

These tests need neither rdkit nor Open Babel: ``rdkit.Chem`` is imported lazily and the
subprocess is stubbed, so what is asserted is the plumbing -- that a configured path reaches the
command line and the error message, and that the default is unchanged for everyone who has
``obabel`` on ``PATH``.
"""

import unittest
from unittest import mock

from pestifer.core.config import Config
from pestifer.charmmff.ligand_paramgen import mol2_writer as MW


class _FakeCommand:
    """Stands in for core.command.Command, recording the command line it was handed."""
    last = None

    def __init__(self, cmdline):
        _FakeCommand.last = cmdline
        self.stderr = 'obabel: command not found'
        self.stdout = ''

    def run(self, **kwargs):
        return 127                      # what a missing binary returns through /bin/sh


class TestObabelIsConfigurable(unittest.TestCase):

    def test_paths_obabel_defaults_to_the_bare_command(self):
        c = Config(userdict={}, quiet=True).configure()
        self.assertEqual(c.shell_commands['obabel'], 'obabel')

    def test_paths_obabel_can_be_relocated(self):
        c = Config(userdict={'paths': {'obabel': '/opt/openbabel/bin/obabel-3.1'}},
                   quiet=True).configure()
        self.assertEqual(c.shell_commands['obabel'], '/opt/openbabel/bin/obabel-3.1')

    def test_obabel_is_not_required_to_load_a_configuration(self):
        """A missing Open Babel must never block a build.

        It is used by one subcommand and by no build task, so requiring it would break pestifer
        for everyone who never runs `make-ligand-mol2`.
        """
        with mock.patch('pestifer.core.config.shutil.which',
                        side_effect=lambda cmd: None if 'obabel' in cmd else f'/usr/bin/{cmd}'):
            c = Config(userdict={}, quiet=True).configure()
        self.assertEqual(c.shell_commands['obabel'], 'obabel',
                         'an unresolvable obabel should be recorded, not fatal')


class TestWriterUsesTheConfiguredPath(unittest.TestCase):

    @staticmethod
    def _mol():
        """An RDKit Mol stand-in with no atoms -- enough to reach the obabel call."""
        m = mock.Mock()
        m.GetAtoms.return_value = []
        m.GetConformer.return_value = mock.Mock()
        return m

    def _write(self, obabel):
        # rdkit is an optional dependency and is imported lazily, so the stand-in only needs the
        # two entry points write_mol2 actually calls before it reaches obabel.
        fake_chem = mock.Mock()
        fake_chem.MolToPDBBlock.return_value = 'ATOM      1  C   LIG A   1\nEND\n'
        with mock.patch.object(MW, 'Command', _FakeCommand), \
             mock.patch.object(MW, '_import_rdkit_chem', return_value=fake_chem):
            with self.assertRaises(MW.Mol2WriteError) as cm:
                MW.write_mol2(self._mol(), 'LIG', '/tmp/does-not-matter.mol2', obabel=obabel)
        return _FakeCommand.last, str(cm.exception)

    def test_default_invocation_is_unchanged(self):
        cmdline, _msg = self._write('obabel')
        self.assertTrue(cmdline.startswith('obabel '), cmdline)

    def test_a_configured_path_reaches_the_command_line(self):
        cmdline, _msg = self._write('/opt/openbabel/bin/obabel-3.1')
        self.assertTrue(cmdline.startswith('/opt/openbabel/bin/obabel-3.1 '), cmdline)

    def test_the_error_names_the_binary_that_was_actually_tried(self):
        # a user who relocated obabel and got the path wrong must be told which path failed,
        # not the generic name
        _cmdline, msg = self._write('/opt/openbabel/bin/obabel-3.1')
        self.assertIn('/opt/openbabel/bin/obabel-3.1', msg)
