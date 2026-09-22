# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Which NAMD a nested build launches.

A build that meets a missing PDB-repository entry runs a whole pipeline of its own inside the same
process, and that pipeline builds its ``Config`` from defaults.  So it resolved ``namd3`` and
``charmrun`` off ``PATH`` while the parent used the absolute paths the user configured under
``paths:``.  Where those differ, the nested build runs a *different binary*: reported 2026-09-22
from Picotte, where a multicore-CUDA NAMD pinned in ``~/.bashrc`` served a CPU-only conformer build
and died with "no CUDA-capable device is detected".

Pure logic over the config: no toolchain needed, and this is exactly the code that misbehaves on a
cluster CI cannot reach.
"""
import os
import unittest
from unittest import mock

from pestifer.core.config import Config


class _Cfg(Config):
    """A Config with just enough state for ``_set_shell_commands``."""

    def __init__(self, paths, namd=None):
        self.data = {'user': {'paths': paths,
                              'namd': {'processor-type': 'cpu', 'deprecated3': {},
                                       **(namd or {})}}}
        self.shell_commands = {}
        self.slurmvars = {}

    def __getitem__(self, key):
        return self.data[key]

    def _verify_catdcd_version(self):
        pass

    def _report_gpu_mode(self):
        pass


DEFAULT_PATHS = {'charmrun': 'charmrun', 'namd3': 'namd3', 'namd2': 'namd2', 'vmd': 'vmd',
                 'catdcd': 'catdcd', 'namd3gpu': 'namd3', 'obabel': 'obabel'}

WHICH = {'charmrun': '/opt/bin/charmrun', 'namd3': '/opt/bin/namd3', 'vmd': '/opt/bin/vmd',
         'catdcd': '/opt/bin/catdcd', 'obabel': '/opt/bin/obabel'}


def _which(cmd):
    """PATH resolution: absolute paths resolve to themselves, bare names to /opt/bin."""
    if cmd.startswith('/'):
        return cmd
    return WHICH.get(cmd)


class TestToolchainInheritance(unittest.TestCase):

    def setUp(self):
        self._saved = Config._resolved_toolchain
        Config._resolved_toolchain = {}

    def tearDown(self):
        Config._resolved_toolchain = self._saved

    def _resolve(self, paths, verify_access=True, processor_type='cpu'):
        c = _Cfg(dict(paths), namd={'processor-type': processor_type})
        with mock.patch('shutil.which', _which), \
             mock.patch.object(os, 'access', return_value=True):
            c._set_shell_commands(verify_access=verify_access)
        return c

    def test_a_configured_path_is_published_for_nested_builds(self):
        parent = self._resolve({**DEFAULT_PATHS, 'namd3': '/ifs/group/namd3'})
        self.assertEqual(parent.shell_commands['namd3'], '/ifs/group/namd3')
        self.assertEqual(Config._resolved_toolchain['namd3'], '/ifs/group/namd3')

    def test_a_nested_build_on_defaults_adopts_it(self):
        """The reported case: the nested build's own paths are defaults, and PATH leads elsewhere."""
        self._resolve({**DEFAULT_PATHS, 'namd3': '/ifs/group/namd3'})
        nested = self._resolve(DEFAULT_PATHS, verify_access=False)
        self.assertEqual(nested.shell_commands['namd3'], '/ifs/group/namd3')

    def test_a_nested_build_that_names_its_own_binary_keeps_it(self):
        self._resolve({**DEFAULT_PATHS, 'namd3': '/ifs/group/namd3'})
        nested = self._resolve({**DEFAULT_PATHS, 'namd3': '/other/namd3'}, verify_access=False)
        self.assertEqual(nested.shell_commands['namd3'], '/other/namd3')

    def test_a_bare_name_is_resolved_to_the_path_it_will_launch(self):
        """The banner printed the absolute path while the launch used the bare name, which mpirun
        re-resolves on the compute node -- so the banner was not a statement about what would run."""
        c = self._resolve(DEFAULT_PATHS)
        self.assertEqual(c.shell_commands['namd3'], '/opt/bin/namd3')
        self.assertEqual(c.shell_commands['charmrun'], '/opt/bin/charmrun')

    def test_an_inherited_namd3_does_not_turn_on_gpu_mode(self):
        """The regression the integration gate caught before 3.23.0 shipped.

        `auto` decides there is a separate GPU binary when `paths.namd3gpu != paths.namd3`, and
        both default to 'namd3'.  Comparing a RESOLVED namd3 against that default made them
        differ, so a CPU box resolved to GPU mode and NAMD was launched with +devices, which it
        rejects outright ("Unknown command-line option +devices").
        """
        self._resolve({**DEFAULT_PATHS, 'namd3': '/ifs/group/namd3'}, processor_type='auto')
        nested = self._resolve(DEFAULT_PATHS, verify_access=False, processor_type='auto')
        self.assertEqual(nested.namd_type, 'cpu')
        self.assertEqual(nested.shell_commands['namd3gpu'], nested.shell_commands['namd3'])

    def test_a_plain_config_on_defaults_is_still_cpu(self):
        c = self._resolve(DEFAULT_PATHS, processor_type='auto')
        self.assertEqual(c.namd_type, 'cpu')

    def test_a_real_separate_gpu_binary_is_still_honoured(self):
        WHICH['namd3gpu'] = '/opt/bin/namd3gpu'
        try:
            c = self._resolve({**DEFAULT_PATHS, 'namd3gpu': 'namd3gpu'}, processor_type='auto')
            self.assertEqual(c.namd_type, 'gpu')
            self.assertEqual(c.shell_commands['namd3gpu'], '/opt/bin/namd3gpu')
        finally:
            WHICH.pop('namd3gpu', None)

    def test_a_post_processing_config_does_not_publish_a_toolchain(self):
        # verify_access=False is a standalone subcommand, not a build; it must not overwrite what a
        # build resolved
        self._resolve({**DEFAULT_PATHS, 'namd3': '/ifs/group/namd3'})
        self._resolve({**DEFAULT_PATHS, 'namd3': '/other/namd3'}, verify_access=False)
        self.assertEqual(Config._resolved_toolchain['namd3'], '/ifs/group/namd3')
