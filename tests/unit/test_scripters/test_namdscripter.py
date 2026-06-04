import unittest
import os
from pestifer.core.config import Config
from pestifer.scripters import NAMDScripter

class TestNAMDScripter(unittest.TestCase):

    def test_charmm(self):
        c = Config().configure_new()
        p: NAMDScripter = c.get_scripter('namd')
        self.assertIsInstance(p, NAMDScripter)
        self.assertEqual(p.namd_version, 3)
        c.RM.charmmff_content.clean_local_charmmff_files()


class TestNAMDLaunchCommand(unittest.TestCase):
    """
    Exercise NAMDScripter._build_launch_command for each CPU launcher mode without
    needing a full Config/charmmff setup. The object is created with __new__ and only
    the attributes that _build_launch_command reads are populated.
    """

    def _make_scripter(self, *, slurmvars, launcher='auto', ncpus=192,
                       namd_type='cpu'):
        p = NAMDScripter.__new__(NAMDScripter)
        p.scriptname = 'job.namd'
        p.namd = 'namd3'
        p.namdgpu = 'namd3'
        p.charmrun = 'charmrun'
        p.namd_type = namd_type
        p.ncpus = ncpus
        p.local_ncpus = 48
        p.ngpus = 0
        p.gpu_devices = ''
        p.slurmvars = slurmvars
        p.namd_config = {'cpu-parallel-launcher': launcher}
        return p

    def test_no_slurm_uses_charmrun(self):
        p = self._make_scripter(slurmvars={})
        c = p._build_launch_command()
        self.assertEqual(c.command, 'charmrun +p 192 namd3 job.namd')
        self.assertTrue(p._single_node_launch)

    def test_auto_single_node_uses_numactl(self):
        p = self._make_scripter(
            slurmvars={'SLURM_NNODES': '1', 'SLURM_NTASKS_PER_NODE': '48'},
            launcher='auto', ncpus=48)
        c = p._build_launch_command()
        self.assertEqual(c.command,
                         'numactl --interleave=all namd3 +p 48 job.namd')
        # single node -> node-local parameter staging is valid
        self.assertTrue(p._single_node_launch)

    def test_auto_multi_node_uses_srun(self):
        p = self._make_scripter(
            slurmvars={'SLURM_NNODES': '4', 'SLURM_NTASKS_PER_NODE': '48'},
            launcher='auto', ncpus=192)
        c = p._build_launch_command()
        self.assertEqual(c.command, 'srun namd3 job.namd')
        # multi-node -> must NOT stage params to node-local scratch
        self.assertFalse(p._single_node_launch)

    def test_explicit_srun(self):
        p = self._make_scripter(
            slurmvars={'SLURM_NNODES': '1', 'SLURM_NTASKS_PER_NODE': '48'},
            launcher='srun')
        c = p._build_launch_command()
        self.assertEqual(c.command, 'srun namd3 job.namd')
        self.assertFalse(p._single_node_launch)

    def test_explicit_srun_with_mpi_type(self):
        p = self._make_scripter(
            slurmvars={'SLURM_NNODES': '4', 'SLURM_NTASKS_PER_NODE': '48'},
            launcher='srun', ncpus=192)
        p.namd_config = {'cpu-parallel-launcher': 'srun', 'srun-mpi-type': 'pmi2'}
        c = p._build_launch_command()
        self.assertEqual(c.command, 'srun --mpi=pmi2 namd3 job.namd')
        self.assertFalse(p._single_node_launch)

    def test_explicit_mpirun(self):
        p = self._make_scripter(
            slurmvars={'SLURM_NNODES': '4', 'SLURM_NTASKS_PER_NODE': '48'},
            launcher='mpirun', ncpus=192)
        c = p._build_launch_command()
        self.assertEqual(c.command, 'mpirun -np 192 namd3 job.namd')
        self.assertFalse(p._single_node_launch)

    def test_explicit_charmrun_uses_mpiexec(self):
        p = self._make_scripter(
            slurmvars={'SLURM_NNODES': '4', 'SLURM_NTASKS_PER_NODE': '48'},
            launcher='charmrun', ncpus=192)
        c = p._build_launch_command()
        self.assertEqual(c.command,
                         'charmrun +p 192 ++mpiexec namd3 job.namd')
        self.assertFalse(p._single_node_launch)

    def test_explicit_numactl_even_multi_node(self):
        p = self._make_scripter(
            slurmvars={'SLURM_NNODES': '4', 'SLURM_NTASKS_PER_NODE': '48'},
            launcher='numactl', ncpus=192)
        c = p._build_launch_command()
        self.assertEqual(c.command,
                         'numactl --interleave=all namd3 +p 192 job.namd')
        # user forced node-local launcher -> staging stays enabled
        self.assertTrue(p._single_node_launch)

    def test_auto_default_when_option_absent(self):
        # namd_config without the key should still default to 'auto'
        p = self._make_scripter(
            slurmvars={'SLURM_NNODES': '4', 'SLURM_NTASKS_PER_NODE': '48'},
            ncpus=192)
        p.namd_config = {}
        c = p._build_launch_command()
        self.assertEqual(c.command, 'srun namd3 job.namd')
        self.assertFalse(p._single_node_launch)
