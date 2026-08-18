# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the build-environment provenance record.

These do not need the external toolchain: the probes are exercised against captured output.
"""

import unittest
from types import SimpleNamespace
from unittest import mock

from pestifer.util import provenance
from pestifer.util.stringthings import __pestifer_version__

_VMD_OUT = """Info) VMD for LINUXAMD64, version 2.0.0 (March 25, 2026)
Info) http://www.ks.uiuc.edu/Research/vmd/
"""

_NAMD_CPU_OUT = """Charm++> Running in Multicore mode: 1 threads (PEs)
Info: NAMD 3.0.2 for Linux-x86_64-multicore
Info: Built Wed Aug 27 15:39:24 CDT 2025 by dhardy on athine.ks.uiuc.edu
ERROR: 'cutoff' is a required configuration option
"""

_NAMD_GPU_OUT = """Info: Built with CUDA version 11080
Info: NAMD 3.0.2 for Linux-x86_64-multicore-CUDA
Info: Built Wed Aug 27 15:48:11 CDT 2025 by dhardy on athine.ks.uiuc.edu
"""


class TestVersionProbes(unittest.TestCase):

    def test_vmd_version_parsed(self):
        with mock.patch.object(provenance, '_run', return_value=_VMD_OUT):
            self.assertEqual(provenance.vmd_version('vmd'), '2.0.0 (March 25, 2026)')

    def test_namd_version_parsed_despite_config_error(self):
        # NAMD prints its banner and *then* rejects the empty probe config; a nonzero exit and
        # an ERROR line must not stop us reading the version.
        with mock.patch.object(provenance, '_run', return_value=_NAMD_CPU_OUT):
            self.assertEqual(provenance.namd_version('namd3'),
                             'NAMD 3.0.2 for Linux-x86_64-multicore')

    def test_namd_version_records_cuda_build(self):
        with mock.patch.object(provenance, '_run', return_value=_NAMD_GPU_OUT):
            self.assertEqual(provenance.namd_version('namd3gpu'),
                             'NAMD 3.0.2 for Linux-x86_64-multicore-CUDA (CUDA 11080)')

    def test_probes_degrade_to_unknown(self):
        # a probe that returns nothing (missing binary, timeout) must not raise
        with mock.patch.object(provenance, '_run', return_value=''):
            self.assertEqual(provenance.vmd_version('vmd'), 'unknown')
            self.assertEqual(provenance.namd_version('namd3'), 'unknown')

    def test_run_never_raises(self):
        self.assertEqual(provenance._run(['definitely-not-a-real-binary-xyzzy']), '')


class TestPackageVersions(unittest.TestCase):

    def test_pestifer_version_comes_from_source_tree_not_metadata(self):
        # importlib.metadata records the version at install time, so an editable install goes
        # stale the moment a release is cut; the record must not inherit that bug.
        with mock.patch('pestifer.util.provenance.importlib.metadata.version',
                        return_value='0.0.0-stale'):
            self.assertEqual(provenance.package_versions()['pestifer'], __pestifer_version__)

    def test_declared_dependencies_are_reported(self):
        pkgs = provenance.package_versions()
        for expected in ('numpy', 'pidibble', 'ycleptic'):
            self.assertIn(expected, pkgs)

    def test_missing_distribution_is_unknown_not_fatal(self):
        with mock.patch.object(provenance, '_dist_names', return_value=['no-such-dist-xyzzy']):
            self.assertEqual(provenance.package_versions()['no-such-dist-xyzzy'], 'unknown')

    def test_extras_are_excluded(self):
        reqs = ['numpy>=1.24', 'pytest; extra == "test"', 'rdkit; extra == "ligand-paramgen"']
        with mock.patch('pestifer.util.provenance.importlib.metadata.requires',
                        return_value=reqs):
            self.assertEqual(provenance._dist_names(), ['numpy'])


class TestExecutableSelection(unittest.TestCase):

    @staticmethod
    def _config(namd_type):
        return SimpleNamespace(namd_type=namd_type,
                               shell_commands={'namd3': 'namd3', 'namd3gpu': 'namd3gpu',
                                               'vmd': 'vmd', 'charmrun': 'charmrun'})

    def _probe(self, namd_type):
        with mock.patch.object(provenance.shutil, 'which', side_effect=lambda c: f'/usr/bin/{c}'), \
             mock.patch.dict(provenance._PROBERS, {}, clear=True):
            return provenance.executable_versions(self._config(namd_type))

    def test_cpu_run_does_not_probe_the_gpu_binary(self):
        # probing the CUDA build binds a GPU device; a CPU run will never invoke it
        self.assertNotIn('namd3gpu', self._probe('cpu'))

    def test_gpu_run_probes_both_binaries(self):
        # a GPU build still runs the CPU binary for any task carrying cpu-override
        got = self._probe('gpu')
        self.assertIn('namd3gpu', got)
        self.assertIn('namd3', got)

    def test_unresolvable_command_recorded_as_unknown(self):
        with mock.patch.object(provenance.shutil, 'which', return_value=None), \
             mock.patch.dict(provenance._PROBERS, {}, clear=True):
            got = provenance.executable_versions(self._config('cpu'))
        self.assertEqual(got['vmd']['path'], 'unknown')


class TestLogEnvironment(unittest.TestCase):

    def test_emits_a_greppable_block(self):
        cfg = SimpleNamespace(namd_type='cpu', shell_commands={})
        lines = []
        with mock.patch.object(provenance, 'charmmff_release', return_value='feb26'):
            provenance.log_environment(cfg, log=lines.append)
        self.assertTrue(all(l.startswith('environment: ') for l in lines))
        joined = '\n'.join(lines)
        self.assertIn('feb26', joined)
        self.assertIn(__pestifer_version__, joined)

    def test_never_raises_when_the_report_fails(self):
        # provenance is worth a second of startup, never a failed build
        with mock.patch.object(provenance, 'environment_report', side_effect=RuntimeError('boom')):
            self.assertEqual(provenance.log_environment(SimpleNamespace()), {})

    def test_charmmff_release_degrades_to_unknown(self):
        self.assertEqual(provenance.charmmff_release(SimpleNamespace()), 'unknown')


class TestProbeFailureIsVisible(unittest.TestCase):
    """An 'unknown' in the record must be explained, not silent: a provenance record that
    degrades without saying so is the failure mode this module exists to prevent."""

    def test_timeout_warns(self):
        import subprocess
        with mock.patch.object(provenance.subprocess, 'run',
                               side_effect=subprocess.TimeoutExpired('vmd', 1)):
            with self.assertLogs('pestifer.util.provenance', level='WARNING') as cm:
                self.assertEqual(provenance._run(['vmd', '--version']), '')
        self.assertIn('timed out', ''.join(cm.output))

    def test_other_failure_warns(self):
        with mock.patch.object(provenance.subprocess, 'run', side_effect=OSError('boom')):
            with self.assertLogs('pestifer.util.provenance', level='WARNING') as cm:
                self.assertEqual(provenance._run(['vmd', '--version']), '')
        self.assertIn('recorded as unknown', ''.join(cm.output))

    def test_timeout_is_generous_enough_for_a_vmd_startup(self):
        # vmd --version performs a full startup; a tight bound silently degrades the record
        self.assertGreaterEqual(provenance._PROBE_TIMEOUT, 60)
