import os
import tempfile
import unittest
from unittest.mock import patch

from pestifer.charmmff.charmmffprm import CharmmParamFile
from pestifer.core.config import Config
from pestifer.core.errors import PestiferBuildError
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
                       namd_type='cpu', namd_config=None):
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
        p.namd_config.update(namd_config or {})
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


class TestVacuumRunPECount(unittest.TestCase):
    """A stage with no periodic cell is clamped to one node's cores, and must say so.

    Keeping the clamp is right -- these stages are small and short, and more ranks would cost more
    in communication than they gain -- but on a 4-node allocation it leaves three nodes idle and
    contradicts an explicit --ncpus.  One membrane build ran 2 of its 215 launches at 48 of 192
    PEs with only a debug line to show for it, and explaining that took an audit of the whole log.
    """

    def _scripter(self, launcher='mpirun', **kw):
        # mpirun, not srun: srun takes its rank count from the SLURM allocation and never sees
        # the PE count, so the clamp cannot reach it (see test_srun_is_not_clamped_at_all)
        return TestNAMDLaunchCommand._make_scripter(
            TestNAMDLaunchCommand(), slurmvars={'SLURM_NNODES': '4',
                                                'SLURM_NTASKS_PER_NODE': '48'},
            launcher=launcher, **kw)

    def test_the_clamp_is_reported_when_it_bites(self):
        p = self._scripter()
        with self.assertLogs('pestifer.scripters.namd', level='INFO') as cm:
            p._build_launch_command(local_execution_only=True)
        out = ''.join(cm.output)
        self.assertIn('48 of 192 PEs', out)
        self.assertIn('vacuum-runs-single-node', out)

    def test_the_clamp_bites_on_a_launcher_that_takes_a_pe_count(self):
        p = self._scripter()
        self.assertEqual(p._build_launch_command(local_execution_only=True).command,
                         'mpirun -np 48 namd3 job.namd')

    def test_the_clamp_can_be_turned_off(self):
        p = self._scripter(namd_config={'vacuum-runs-single-node': False})
        with self.assertNoLogs('pestifer.scripters.namd', level='INFO'):
            c = p._build_launch_command(local_execution_only=True)
        self.assertEqual(c.command, 'mpirun -np 192 namd3 job.namd')

    def test_srun_is_not_clamped_at_all(self):
        """srun spawns one rank per allocated task, so the PE count never reaches it -- claiming
        a clamp there would report something that did not happen."""
        p = self._scripter(launcher='srun')
        with self.assertNoLogs('pestifer.scripters.namd', level='INFO'):
            c = p._build_launch_command(local_execution_only=True)
        self.assertEqual(c.command, 'srun namd3 job.namd')

    def test_a_single_node_allocation_is_silent(self):
        # nothing is being clamped away there, so there is nothing to explain
        p = self._scripter(ncpus=48)
        with self.assertNoLogs('pestifer.scripters.namd', level='INFO'):
            p._build_launch_command(local_execution_only=True)

    def test_a_periodic_stage_uses_the_whole_allocation(self):
        p = self._scripter()
        with self.assertNoLogs('pestifer.scripters.namd', level='INFO'):
            c = p._build_launch_command()
        self.assertEqual(c.command, 'mpirun -np 192 namd3 job.namd')


class TestConsolidateParams(unittest.TestCase):
    """``consolidate_params`` reduces the parameter set to one minimal .prm and points the
    script's single ``parameters`` line at it.

    When the set is *already* one consolidated file, the file keeps its own name: packaging
    writes the package's NAMD config under the package basename while shipping the minimal
    .prm generated under the terminate basename, so re-deriving the name from the basename in
    force left the config referencing a file the tarball did not contain.
    """

    PSF = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'inputs',
                                       'existing.psf'))

    def _make_scripter(self, tmp, basename, parameters):
        p = NAMDScripter.__new__(NAMDScripter)
        p.basename = basename
        p.scriptname = os.path.join(tmp, f'{basename}.namd')
        p.parameters = parameters
        with open(p.scriptname, 'w') as fh:
            fh.write('structure sys.psf\n')
            for q in parameters:
                fh.write(f'parameters {q}\n')
            fh.write('minimize 100\n')
        return p

    @staticmethod
    def _param_lines(scriptname):
        return [l.split()[1] for l in open(scriptname).read().splitlines()
                if l.startswith('parameters ')]

    def test_already_minimal_keeps_its_own_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            prm = os.path.join(tmp, 'mysys_minimal.prm')
            open(prm, 'w').write('* minimal\n*\n\nEND\n')
            # scripter basename is the *package* basename; the parameter file is not
            p = self._make_scripter(tmp, 'prod_mysys', [prm])
            out = p.consolidate_params(self.PSF)
            self.assertEqual(out, prm)
            self.assertEqual(self._param_lines(p.scriptname), [prm])
            self.assertFalse(os.path.exists(os.path.join(tmp, 'prod_mysys_minimal.prm')))

    def test_returns_none_without_parameters_or_psf(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = self._make_scripter(tmp, 'mysys', [])
            self.assertIsNone(p.consolidate_params(self.PSF))
            p = self._make_scripter(tmp, 'mysys', ['some.prm'])
            self.assertIsNone(p.consolidate_params(os.path.join(tmp, 'nope.psf')))


class TestConsolidateParamsCollision(unittest.TestCase):
    """Sibling sub-builds restart task numbering, so two of them name this file identically and
    the last writer wins.  Downstream that is either a loud NAMD failure on a missing vdW
    parameter, or -- where the surviving file is a superset drawn from a different merged source
    set -- no failure at all and quietly different values.  Warn where both names are in hand."""

    PSF = TestConsolidateParams.PSF

    def _prm(self, path, types):
        with open(path, 'w') as fh:
            fh.write('* test\n*\n\nNONBONDED\n')
            for t in types:
                fh.write(f'{t}  0.0  -0.1  1.9\n')
            fh.write('\nEND\n')
        return path

    def test_overwriting_a_differing_file_warns(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = self._prm(os.path.join(tmp, 'src.prm'),
                            ['CT1', 'CT2', 'CT3', 'HA1', 'HA2', 'HA3', 'NH1', 'O', 'C', 'OT'])
            p = TestConsolidateParams()._make_scripter(tmp, 'shared', [src])
            cwd = os.getcwd()
            try:
                os.chdir(tmp)
                # a previous sub-build already wrote this name, with a different content
                self._prm('shared_minimal.prm', ['CT1'])
                with self.assertLogs('pestifer.scripters.namd', level='WARNING') as cm:
                    # this .prm carries vdW records only, so it genuinely cannot cover the PSF's
                    # bonded terms; the overwrite warning is emitted before that check fires
                    with self.assertRaises(PestiferBuildError):
                        p.consolidate_params(self.PSF)
                self.assertTrue(any('sharing one artifact name' in m for m in cm.output),
                                cm.output)
            finally:
                os.chdir(cwd)


# A complete parameter set for the 4-atom PSF below: every bond, angle, dihedral, improper
# and vdW record its topology needs, and nothing else.
_COMPLETE_PRM = """\
* complete synthetic parameters
*

BONDS
NH1  CT1   300.0   1.45
CT1  C     250.0   1.49
C    O     620.0   1.23

ANGLES
NH1  CT1  C    50.0   110.0
CT1  C    O    80.0   121.0

DIHEDRALS
X    CT1  C    X     0.2000  1   0.00

IMPROPER
C    X    X    O    120.0   0   0.00

NONBONDED nbxmod 5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 16.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
NH1   0.0  -0.20   1.85
CT1   0.0  -0.032  2.00
C     0.0  -0.11   2.00
O     0.0  -0.12   1.70
"""

_TINY_PSF = """\
PSF EXT

       3 !NTITLE
 REMARKS tiny test system
 REMARKS topology top_all36_prot.rtf
 REMARKS segment PROA { first NTER; last CTER; auto angles dihedrals }

       4 !NATOM
       1 PROA     7        THR       N        NH1   -0.470000       14.0070           0
       2 PROA     7        THR       CA       CT1    0.070000       12.0110           0
       3 PROA     7        THR       C        C      0.510000       12.0110           0
       4 PROA     7        THR       O        O     -0.510000       15.9990           0

       3 !NBOND: bonds
       1       2       2       3       3       4

       2 !NTHETA: angles
       1       2       3       2       3       4

       1 !NPHI: dihedrals
       1       2       3       4

       1 !NIMPHI: impropers
       3       1       2       4

       0 !NDON: donors

       0 !NACC: acceptors

       0 !NNB
"""


class TestConsolidateParamsVerifiesCoverage(unittest.TestCase):
    """A parameter the PSF needs and the run does not have is fatal to NAMD either way.

    Checking it where the minimal .prm is written turns it into an error that names the residue
    and the term, instead of a NAMD abort that names only atom serials.  Validated against 63
    real builds: 62 clean, and the one flagged (a THR whose CB was retyped CT2 by a patch that
    left HB as HA1) is a build NAMD really did kill with
    ``UNABLE TO FIND ANGLE PARAMETERS FOR CT1 CT2 HA1 (ATOMS 3790 3792 3794)``.
    """

    def _setup(self, tmp, prm_text):
        psf = os.path.join(tmp, 'tiny.psf')
        open(psf, 'w').write(_TINY_PSF)
        prm = os.path.join(tmp, 'src.prm')
        open(prm, 'w').write(prm_text)
        p = NAMDScripter.__new__(NAMDScripter)
        p.basename = 'tiny'
        p.scriptname = os.path.join(tmp, 'tiny.namd')
        p.parameters = [prm]
        with open(p.scriptname, 'w') as fh:
            fh.write(f'structure tiny.psf\nparameters {prm}\nminimize 100\n')
        return p, psf

    def test_complete_set_consolidates_without_complaint(self):
        with tempfile.TemporaryDirectory() as tmp:
            p, psf = self._setup(tmp, _COMPLETE_PRM)
            cwd = os.getcwd()
            try:
                os.chdir(tmp)
                out = p.consolidate_params(psf)
                self.assertEqual(out, 'tiny_minimal.prm')
                self.assertTrue(os.path.exists(out))
            finally:
                os.chdir(cwd)

    def test_a_term_absent_from_every_source_file_raises_naming_the_residue(self):
        """The failure the check exists for, exercised through consolidate_params itself --
        not through the helper -- so removing the call makes this test go red."""
        with tempfile.TemporaryDirectory() as tmp:
            p, psf = self._setup(tmp, _COMPLETE_PRM.replace('CT1  C    O    80.0   121.0', ''))
            cwd = os.getcwd()
            try:
                os.chdir(tmp)
                with self.assertRaises(PestiferBuildError) as cm:
                    p.consolidate_params(psf)
            finally:
                os.chdir(cwd)
            msg = str(cm.exception)
            self.assertIn('CT1-C-O', msg)
            self.assertIn('THR PROA7', msg)      # names the residue NAMD would not
            self.assertIn('charmmff.standard.str', msg)
            self.assertNotIn('PESTIFER BUG', msg)

    def test_a_term_lost_in_consolidation_is_reported_as_a_pestifer_bug(self):
        """Present in the merged set but absent from the consolidated file is the opposite
        failure -- an extraction bug -- and must not be reported as a user misconfiguration."""
        real = CharmmParamFile.extract_for_atomtypes

        def lossy(self, atomtypes):
            out = real(self, atomtypes)
            out.angles = [a for a in out.angles
                          if (a.type1, a.type2, a.type3) != ('CT1', 'C', 'O')]
            return out

        with tempfile.TemporaryDirectory() as tmp:
            p, psf = self._setup(tmp, _COMPLETE_PRM)
            cwd = os.getcwd()
            try:
                os.chdir(tmp)
                with patch.object(CharmmParamFile, 'extract_for_atomtypes', lossy):
                    with self.assertRaises(PestiferBuildError) as cm:
                        p.consolidate_params(psf)
            finally:
                os.chdir(cwd)
            msg = str(cm.exception)
            self.assertIn('PESTIFER BUG', msg)
            self.assertIn('CT1-C-O', msg)
            self.assertNotIn('charmmff.standard.str', msg)
