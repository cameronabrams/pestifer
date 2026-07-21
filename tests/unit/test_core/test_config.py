# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import pytest
# these tests construct a verify_access Config, which requires the external
# vmd/namd/charmrun toolchain; skip them where those binaries are absent
pytestmark = pytest.mark.needs_tools

import unittest
import os
import shutil
import subprocess
from unittest import mock
import yaml
from pestifer.core.config import Config, detect_local_gpu_ids
from pestifer.core.labels import Labels

class TestConfig(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.c = Config(quiet=True).configure_new()
        cls.RM = cls.c.RM

    def test_labels(self):
        c = self.c
        self.assertEqual(os.path.join(c.tcl_root,'vmdrc.tcl'), c.vmd_startup_script)
        self.assertEqual(Labels.segtype_of_resname['ALA'], 'protein')
        self.assertEqual(Labels.segtype_of_resname['MAN'], 'glycan')
        self.assertEqual(Labels.segtype_of_resname['CL'], 'ion')
        self.assertEqual(Labels.charmm_resname_of_pdb_resname.get('MAN'), 'AMAN')
        self.assertEqual(Labels.charmm_resname_of_pdb_resname.get('ALA'), None)
        self.assertEqual(Labels.res_123['A'], 'ALA')
        self.assertEqual(Labels.res_321['PHE'], 'F')

    def test_config_userdict(self):
        C = self.c
        ud = C['user']
        D = Config(userdict=ud).configure_new()
        self.assertTrue('user' in D)
        
    def test_detect_local_gpu_ids_no_nvidia_smi(self):
        # nvidia-smi not on PATH -> empty list (CPU-only workstation)
        with mock.patch('pestifer.core.config.shutil.which', return_value=None):
            self.assertEqual(detect_local_gpu_ids(), [])

    def test_detect_local_gpu_ids_parses_indices(self):
        fake = subprocess.CompletedProcess(args=[], returncode=0,
                                           stdout='0\n1\n2\n3\n', stderr='')
        with mock.patch('pestifer.core.config.shutil.which', return_value='/usr/bin/nvidia-smi'), \
             mock.patch('pestifer.core.config.subprocess.run', return_value=fake):
            self.assertEqual(detect_local_gpu_ids(), [0, 1, 2, 3])

    def test_detect_local_gpu_ids_handles_failure(self):
        # nvidia-smi present but errors (e.g. driver/library mismatch) -> empty list
        fake = subprocess.CompletedProcess(args=[], returncode=9,
                                           stdout='', stderr='boom')
        with mock.patch('pestifer.core.config.shutil.which', return_value='/usr/bin/nvidia-smi'), \
             mock.patch('pestifer.core.config.subprocess.run', return_value=fake):
            self.assertEqual(detect_local_gpu_ids(), [])
        # also tolerate the executable raising
        with mock.patch('pestifer.core.config.shutil.which', return_value='/usr/bin/nvidia-smi'), \
             mock.patch('pestifer.core.config.subprocess.run',
                        side_effect=OSError('nope')):
            self.assertEqual(detect_local_gpu_ids(), [])

    # --- GPU-not-elected warning -------------------------------------------------
    # pestifer enables GPU mode only when paths.namd3gpu differs from paths.namd3; both
    # default to 'namd3', so on a host whose GPU-resident build is a separate binary the
    # GPU is silently unused.  These pin the warn/no-warn logic without needing a real GPU.

    @staticmethod
    def _gpu_stub(ngpus, quiet=False):
        """A bare Config for exercising _warn_gpu_not_elected in isolation."""
        c = Config.__new__(Config)
        c.ngpus = ngpus
        c.quiet = quiet
        return c

    @staticmethod
    def _which(mapping):
        return lambda name: mapping.get(name)

    def test_warn_gpu_not_elected_fires(self):
        # GPU present + separate namd3gpu on PATH + paths.namd3gpu still == paths.namd3
        c = self._gpu_stub(ngpus=1)
        which = self._which({'namd3gpu': '/usr/local/bin/namd3gpu',
                             'namd3': '/usr/local/bin/namd3'})
        with mock.patch('pestifer.core.config.shutil.which', side_effect=which):
            with self.assertLogs('pestifer.core.config', level='WARNING') as cm:
                c._warn_gpu_not_elected('namd3', 'namd3')
        msg = ''.join(cm.output)
        self.assertIn('/usr/local/bin/namd3gpu', msg)
        self.assertIn('GPU mode is OFF', msg)

    def test_no_warn_when_no_gpu(self):
        c = self._gpu_stub(ngpus=0)
        which = self._which({'namd3gpu': '/usr/local/bin/namd3gpu',
                             'namd3': '/usr/local/bin/namd3'})
        with mock.patch('pestifer.core.config.shutil.which', side_effect=which):
            with self.assertNoLogs('pestifer.core.config', level='WARNING'):
                c._warn_gpu_not_elected('namd3', 'namd3')

    def test_no_warn_when_no_separate_gpu_build(self):
        # GPU present but no GPU-resident binary to point at
        c = self._gpu_stub(ngpus=1)
        which = self._which({'namd3': '/usr/local/bin/namd3'})
        with mock.patch('pestifer.core.config.shutil.which', side_effect=which):
            with self.assertNoLogs('pestifer.core.config', level='WARNING'):
                c._warn_gpu_not_elected('namd3', 'namd3')

    def test_no_warn_when_single_binary_handles_both(self):
        # namd3gpu resolves to the same binary as namd3 -- nothing to elect
        c = self._gpu_stub(ngpus=1)
        which = self._which({'namd3gpu': '/usr/local/bin/namd3',
                             'namd3': '/usr/local/bin/namd3'})
        with mock.patch('pestifer.core.config.shutil.which', side_effect=which):
            with self.assertNoLogs('pestifer.core.config', level='WARNING'):
                c._warn_gpu_not_elected('namd3', 'namd3')

    def test_no_warn_when_quiet(self):
        # subcontroller configs (taskless_subconfig) are quiet; don't duplicate the warning
        c = self._gpu_stub(ngpus=1, quiet=True)
        which = self._which({'namd3gpu': '/usr/local/bin/namd3gpu',
                             'namd3': '/usr/local/bin/namd3'})
        with mock.patch('pestifer.core.config.shutil.which', side_effect=which):
            with self.assertNoLogs('pestifer.core.config', level='WARNING'):
                c._warn_gpu_not_elected('namd3', 'namd3')

    def test_taskless_subconfig_inherits_user_namd(self):
        # A subcontroller config (e.g. for make_membrane_system relaxation MD) must
        # inherit the parent's NAMD settings, not fall back to schema defaults.
        parent = Config(userdict={'namd': {'cpu-parallel-launcher': 'mpirun'}},
                        quiet=True, RM=self.RM).configure()
        self.assertEqual(parent['user']['namd']['cpu-parallel-launcher'], 'mpirun')
        sub = parent.taskless_subconfig()
        self.assertEqual(sub['user']['namd']['cpu-parallel-launcher'], 'mpirun')
        self.assertEqual(len(sub['user']['tasks']), 0)
        # the NAMD scripter the subcontroller would use must see it too
        self.assertEqual(sub.get_scripter('namd').namd_config['cpu-parallel-launcher'],
                         'mpirun')

    def test_config_user(self):
        tmpdir = '__test_config_user'
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        os.chdir(tmpdir)
        EM = self.RM.example_manager
        e1 = EM.checkout_example(1)
        c = Config(userfile=e1.scriptname).configure_new()
        self.assertTrue('user' in c)
        self.assertTrue('tasks' in c['user'])
        self.assertEqual(len(c['user']['tasks'][0].keys()), 1)
        task1 = c['user']['tasks'][0]
        name, specs = [(x, y) for x, y in task1.items()][0]
        self.assertEqual(name, 'fetch')
        self.assertTrue('sourceID' in specs)
        self.assertEqual(specs['sourceID'], '6pti')
        task2 = c['user']['tasks'][1]
        name, specs = [(x, y) for x, y in task2.items()][0]
        self.assertEqual(name, 'psfgen')
        task3 = c['user']['tasks'][2]
        name, specs = [(x, y) for x, y in task3.items()][0]
        self.assertEqual(name, 'validate')
        os.chdir('..')

    def test_config_task_source(self):
        tmpdir = '__test_config_task_source'
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        os.chdir(tmpdir)
        EM = self.RM.example_manager
        e7 = EM.checkout_example(7)
        c = Config(userfile=e7.scriptname).configure_new()
        self.assertTrue('user' in c)
        self.assertTrue('tasks' in c['user'])
        self.assertEqual(len(c['user']['tasks'][0].keys()), 1)
        task1 = c['user']['tasks'][0]
        name, specs = [(x, y) for x, y in task1.items()][0]
        self.assertEqual(name, 'fetch')
        self.assertEqual(specs['sourceID'], '4zmj')
        task2 = c['user']['tasks'][1]
        name, specs = [(x, y) for x, y in task2.items()][0]
        self.assertEqual(name, 'psfgen')
        sourcespecs = specs['source']
        self.assertTrue('transform_reserves' in sourcespecs)
        resm = sourcespecs['transform_reserves']
        self.assertTrue('G' in resm)
        self.assertTrue('B' in resm)
        self.assertTrue(resm['G'] == ['H', 'J'])
        self.assertTrue(resm['B'] == ['C', 'D'])
        os.chdir('..')

    def test_config_help_examples(self):
        EM = self.RM.example_manager
        tmpdir = '__test_config_help_examples'
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.mkdir(tmpdir)
        os.chdir(tmpdir)
        try:
            # Every example must check out and expose the base+user schema correctly.
            # EM.examples is in filesystem order, not id order, so drive each checkout by
            # the example's own id (never by enumerate position).
            for example in EM.examples:
                EM.checkout_example(example.example_id)
                configfile_path = os.path.join(EM.path, example.scriptpath)
                self.assertTrue(os.path.exists(configfile_path), f'Config file {configfile_path} does not exist after checkout')
                c = Config(userfile=configfile_path).configure_new()
                self.assertTrue('base' in c)
                self.assertTrue('user' in c)
                D = c['base']['attributes']
                self.assertEqual(len(D), 6)
                tld = [x['name'] for x in D]
                self.assertEqual(tld, ['charmmff', 'psfgen', 'namd', 'title', 'paths', 'tasks'])
                T = D[tld.index('title')]
                self.assertEqual(T['name'], 'title')
                self.assertEqual(T['type'], 'str')
                self.assertTrue('text' in T)
                self.assertEqual(T['default'], 'Pestifer')
                U = c['user']
                self.assertTrue('title' in U)
                self.assertTrue('paths' in U)
                self.assertTrue('tasks' in U)
                self.assertIsInstance(U['tasks'], list)
                self.assertTrue(len(U['tasks']) > 0, f'{example.shortname}: no tasks')
                paths = U['paths']
                # An example may legitimately override any of these defaults with an
                # absolute path. Just check the resolved entry ends in the expected
                # binary name.
                self.assertTrue(paths['namd3'].endswith('namd3'))
                self.assertTrue(paths['charmrun'].endswith('charmrun'))
                self.assertTrue(paths['vmd'].endswith('vmd'))
                self.assertTrue(paths['catdcd'].endswith('catdcd'))
                self.assertFalse('charmff' in paths)

            # Not every example fetches or runs psfgen (e.g. the multi-script merge
            # example builds from local inputs), so exercise the fetch + psfgen with
            # mods.ssbondsdelete structure against a canonical example that has it: bpti3.
            bpti3 = next((e for e in EM.examples if e.shortname == 'bpti3'), None)
            self.assertIsNotNone(bpti3, 'expected example bpti3 to be present')
            c = Config(userfile=os.path.join(EM.path, bpti3.scriptpath)).configure_new()
            tasks = c['user']['tasks']
            self.assertTrue(any('fetch' in t for t in tasks))
            psfgen_task = next((t for t in tasks if 'psfgen' in t), None)
            self.assertIsNotNone(psfgen_task, 'No psfgen task found in bpti3')
            specs = psfgen_task['psfgen']
            self.assertTrue('source' in specs)
            self.assertTrue('mods' in specs)
            self.assertTrue('ssbondsdelete' in specs['mods'])
        finally:
            os.chdir('..')
