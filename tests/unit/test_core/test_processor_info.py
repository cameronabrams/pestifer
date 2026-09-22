# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
CPU/GPU resource detection.  Deliberately NOT marked ``needs_tools``: this is pure logic over
environment variables, and it is exactly the code that misbehaves on a cluster nobody can run
CI on, so it must run everywhere.
"""
import os
import unittest
from unittest import mock

from pestifer.core.config import Config


def _cfg(slurmvars):
    c = Config.__new__(Config)          # no toolchain, no config file
    c.slurmvars = slurmvars
    return c


class TestUsableNcpus(unittest.TestCase):
    """
    ``os.cpu_count()`` reports the machine, not the allocation.  Inside a 24-core SLURM job on
    a 48-core node it says 48 and Charm++ then refuses the launch; and computing PEs from
    *tasks* rather than cpus gives 1 for ``--ntasks=1 --cpus-per-task=24``, which runs correctly
    and about 24x slow.
    """

    def test_slurm_cpus_per_task_beats_the_node_count(self):
        with mock.patch.object(os, 'sched_getaffinity', return_value=set(range(48))):
            self.assertEqual(_cfg({'SLURM_CPUS_PER_TASK': '24'})._usable_ncpus(), 24)

    def test_slurm_cpus_on_node_used_when_per_task_absent(self):
        with mock.patch.object(os, 'sched_getaffinity', return_value=set(range(48))):
            self.assertEqual(_cfg({'SLURM_CPUS_ON_NODE': '8'})._usable_ncpus(), 8)

    def test_affinity_wins_when_it_is_the_smaller_bound(self):
        # a site that grants 48 but pins us to 24 -- each source is an upper bound, so take min
        with mock.patch.object(os, 'sched_getaffinity', return_value=set(range(24))):
            self.assertEqual(_cfg({'SLURM_CPUS_PER_TASK': '48'})._usable_ncpus(), 24)

    def test_affinity_used_outside_slurm(self):
        with mock.patch.object(os, 'sched_getaffinity', return_value=set(range(6))):
            self.assertEqual(_cfg({})._usable_ncpus(), 6)

    def test_unparseable_slurm_value_is_ignored_not_fatal(self):
        with mock.patch.object(os, 'sched_getaffinity', return_value=set(range(4))):
            self.assertEqual(_cfg({'SLURM_CPUS_PER_TASK': 'not-a-number'})._usable_ncpus(), 4)

    def test_zero_or_negative_slurm_value_is_ignored(self):
        with mock.patch.object(os, 'sched_getaffinity', return_value=set(range(4))):
            self.assertEqual(_cfg({'SLURM_CPUS_PER_TASK': '0'})._usable_ncpus(), 4)

    def test_falls_back_to_cpu_count_without_sched_getaffinity(self):
        # macOS has no sched_getaffinity
        with mock.patch.object(os, 'cpu_count', return_value=11):
            with mock.patch('pestifer.core.config.os.sched_getaffinity', side_effect=AttributeError):
                self.assertEqual(_cfg({})._usable_ncpus(), 11)

    def test_a_real_call_returns_something_sane(self):
        n = _cfg({})._usable_ncpus()
        self.assertIsInstance(n, int)
        self.assertGreaterEqual(n, 1)


class TestSlurmGpuCount(unittest.TestCase):
    """
    ``ngpus`` was computed from ``gpus_allocated``, which is initialized to '' and never
    assigned -- so ``''.split(',')`` made it 1 for every allocation, however many GPUs.
    """

    def _run_under_slurm(self, job_gpus):
        c = Config.__new__(Config)
        c.data = {'user': {'namd': {'ncpus': 0}}}     # Yclept is a UserDict
        c.ncpus_override = 0
        env = {'SLURM_JOB_ID': '1', 'SLURM_NNODES': '1', 'SLURM_JOB_GPUS': job_gpus,
               'SLURM_CPUS_ON_NODE': '8'}
        # _set_processor_info rebuilds slurmvars from os.environ, so the environment is what
        # has to be faked -- injecting c.slurmvars is silently discarded.
        with mock.patch.dict(os.environ, env, clear=True):
            c._set_processor_info()
        return c.ngpus, c.gpu_devices

    def test_four_gpus_are_counted_as_four(self):
        self.assertEqual(self._run_under_slurm('0,1,2,3'), (4, '0,1,2,3'))

    def test_one_gpu_is_counted_as_one(self):
        self.assertEqual(self._run_under_slurm('0'), (1, '0'))

    def _banner(self, env, affinity=48, cpu_count=48):
        c = Config.__new__(Config)
        c.data = {'user': {'namd': {'ncpus': 0}}}
        c.ncpus_override = 0
        with mock.patch.dict(os.environ, env, clear=True):
            with mock.patch.object(os, 'sched_getaffinity', return_value=set(range(affinity))), \
                 mock.patch.object(os, 'cpu_count', return_value=cpu_count):
                retstr = c._set_processor_info()
        return c, retstr

    def test_a_multi_node_allocation_counts_every_node(self):
        """The regression this exists for: from 3.22.2 the auto-detected count was one node's
        worth, so a 4-node job ran a 25-hour equilibration on 48 of 192 cores with the banner
        reporting "48 cpus" and nothing saying otherwise."""
        c, retstr = self._banner({'SLURM_JOB_ID': '1', 'SLURM_NNODES': '4',
                                  'SLURM_CPUS_ON_NODE': '48',
                                  'SLURM_JOB_CPUS_PER_NODE': '48(x4)'})
        self.assertEqual(c.ncpus, 192, f'got {retstr!r}')
        self.assertEqual(c.local_ncpus, 48, 'the per-node count still feeds +p')
        self.assertIn('192 total', retstr)

    def test_a_heterogeneous_allocation_sums_its_groups(self):
        c, _ = self._banner({'SLURM_JOB_ID': '1', 'SLURM_NNODES': '3',
                             'SLURM_JOB_CPUS_PER_NODE': '72,48(x2)'})
        self.assertEqual(c.ncpus, 168)

    def test_cpus_per_task_does_not_shrink_the_allocation(self):
        """--cpus-per-task=1 made the banner read "1 cpus"; the allocation is still 2x48."""
        c, retstr = self._banner({'SLURM_JOB_ID': '1', 'SLURM_NNODES': '2',
                                  'SLURM_CPUS_PER_TASK': '1', 'SLURM_CPUS_ON_NODE': '48',
                                  'SLURM_JOB_CPUS_PER_NODE': '48(x2)'})
        self.assertEqual(c.ncpus, 96, f'got {retstr!r}')

    def test_an_unparseable_breakdown_falls_back_to_nodes_times_node_count(self):
        c, _ = self._banner({'SLURM_JOB_ID': '1', 'SLURM_NNODES': '2',
                             'SLURM_CPUS_ON_NODE': '48',
                             'SLURM_JOB_CPUS_PER_NODE': 'garbage'})
        self.assertEqual(c.ncpus, 96)

    def test_an_explicit_ncpus_still_wins(self):
        c = Config.__new__(Config)
        c.data = {'user': {'namd': {'ncpus': 0}}}
        c.ncpus_override = 480
        env = {'SLURM_JOB_ID': '1', 'SLURM_NNODES': '10', 'SLURM_JOB_CPUS_PER_NODE': '48(x10)'}
        with mock.patch.dict(os.environ, env, clear=True):
            with mock.patch.object(os, 'sched_getaffinity', return_value=set(range(48))):
                retstr = c._set_processor_info()
        self.assertEqual(c.ncpus, 480)
        self.assertIn('override', retstr)

    def test_cpus_come_from_the_allocation_not_the_node(self):
        c = Config.__new__(Config)
        c.data = {'user': {'namd': {'ncpus': 0}}}
        c.ncpus_override = 0
        env = {'SLURM_JOB_ID': '1', 'SLURM_NNODES': '1', 'SLURM_CPUS_PER_TASK': '24'}
        # cpu_count must be mocked too, and to something != 24: on a 24-core dev box the old
        # os.cpu_count() implementation returns the right answer by coincidence and the test
        # passes against the bug it is meant to catch.
        with mock.patch.dict(os.environ, env, clear=True):
            with mock.patch.object(os, 'sched_getaffinity', return_value=set(range(48))), \
                 mock.patch.object(os, 'cpu_count', return_value=48):
                retstr = c._set_processor_info()
        self.assertEqual(c.ncpus, 24, f'got {retstr!r}')
        self.assertIn('24 cpus', retstr)
