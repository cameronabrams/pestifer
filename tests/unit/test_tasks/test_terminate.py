import os
import shutil
import tarfile
import tempfile
import unittest
from unittest import mock

from pestifer.core.artifacts import (CharmmffParFileArtifact, FileArtifactList, NAMDCoorFileArtifact,
                                     NAMDVelFileArtifact, PDBFileArtifact, PSFFileArtifact,
                                     StateArtifacts)
from pestifer.core.config import Config
from pestifer.core.pipeline import PipelineContext
from pestifer.tasks.terminate import TerminateTask

PSF = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'inputs', 'existing.psf'))

#: NAMD config keywords whose value names a file the packaged tarball must therefore contain.
NAMD_FILE_KEYS = {'structure', 'coordinates', 'bincoordinates', 'binvelocities', 'extendedsystem',
                  'parameters', 'consref', 'conskfile', 'colvarsconfig', 'extrabondsfile'}


class _FA:
    def __init__(self, name):
        self.name = name

    def exists(self):
        return True


class _State:
    def __init__(self, psf):
        self.psf = _FA(psf)
        self.pdb = self.coor = self.xsc = self.vel = None


class TestArchiveName(unittest.TestCase):
    """The build-intermediate tarball defaults to '{basename}-artifacts' so successive builds
    in one directory do not clobber a shared 'artifacts.tar.gz'; an explicit 'artifacts' spec
    still wins."""

    def _task(self, specs):
        t = TerminateTask.__new__(TerminateTask)
        t.specs = specs
        return t

    def test_derives_from_basename(self):
        self.assertEqual(self._task({'basename': 'mysys'})._archive_name(), 'mysys-artifacts')

    def test_explicit_artifacts_wins(self):
        self.assertEqual(self._task({'basename': 'mysys', 'artifacts': 'myarch'})._archive_name(), 'myarch')

    def test_falls_back_without_basename(self):
        self.assertEqual(self._task({})._archive_name(), 'artifacts')
        self.assertEqual(self._task({'artifacts': None, 'basename': None})._archive_name(), 'artifacts')


class TestChainmapFailureIsNotFatal(unittest.TestCase):
    """A 1.58M-atom build lost its package, minimal parameters and run record when the chain map
    could not be written -- the PDB reader cannot read the ``*****`` overflow serials VMD writes
    past 1,048,575 atoms.  Everything in ``do`` after the chain map is independent of it."""

    def _task(self):
        t = TerminateTask.__new__(TerminateTask)
        t.specs = {'chainmapfile': 'cm.yaml'}
        t.next_basename = lambda: None
        t.test_standard = lambda: 0
        t.generate_minimal_params = lambda: None
        t.print_system_report = lambda: None
        t.capture_system_facts = lambda: setattr(t, 'facts_captured', True)
        t.make_package = lambda: 0
        t.cleanup = lambda: 0
        return t

    def test_a_failed_chainmap_does_not_abort_the_task(self):
        t = self._task()

        def boom():
            raise ValueError("invalid literal for int() with base 16: '*****'")

        t.write_chainmaps = boom
        with self.assertLogs('pestifer.tasks.terminate', level='WARNING') as cm:
            result = t.do()
        self.assertEqual(result, 0)
        self.assertTrue(t.facts_captured)          # packaging and the run record still happen
        self.assertIn('*****', ''.join(cm.output))

    def test_a_working_chainmap_is_still_written(self):
        t = self._task()
        t.write_chainmaps = lambda: setattr(t, 'chainmap_written', True)
        self.assertEqual(t.do(), 0)
        self.assertTrue(t.chainmap_written)


class TestResumedRunDoesNotSweep(unittest.TestCase):
    """On a --restart the skipped tasks register nothing in this process's pipeline, so the sweep
    saw 26 files of a 1.58M-atom build: it left ~172 GB loose (the 24 GB production trajectory
    among it) and archived-and-removed the five restored state files instead."""

    def _task(self, resumed_from):
        t = TerminateTask.__new__(TerminateTask)
        t.specs = {}
        t.run_resumed_from = resumed_from

        class _Pipeline:
            def get_all_file_artifacts(self):
                return FileArtifactList([])

        t.pipeline = _Pipeline()
        return t

    def _cleanup(self, t):
        """Run cleanup with the tarball write itself stubbed out; returns the calls it made."""
        swept = []
        with mock.patch.object(FileArtifactList, 'make_tarball',
                               lambda self, name, **kw: swept.append((name, kw.get('remove')))):
            self.assertEqual(t.cleanup(), 0)
        return swept

    def test_a_resumed_run_leaves_the_directory_alone(self):
        t = self._task(13)
        with self.assertLogs('pestifer.tasks.terminate', level='WARNING') as cm:
            swept = self._cleanup(t)
        self.assertEqual(swept, [])
        self.assertIn('skipping the intermediate-file sweep', ''.join(cm.output))

    def test_a_one_pass_run_still_sweeps(self):
        # the sweep removes what it archives, which is right only when this process built it all
        self.assertEqual(self._cleanup(self._task(0)), [('artifacts', True)])


class TestSystemReport(unittest.TestCase):
    def test_reports_total_charge(self):
        t = TerminateTask.__new__(TerminateTask)
        t.get_current_artifact = lambda k: _State(PSF) if k == 'state' else None
        with self.assertLogs('pestifer.tasks.terminate', level='INFO') as cm:
            t.print_system_report()
        out = '\n'.join(cm.output)
        self.assertIn('Total charge', out)
        self.assertIn('5.0000 e', out)          # existing.psf nets to +5 e


class TestPackagedConfigIsSelfContained(unittest.TestCase):
    """Every file the packaged NAMD config names must be in the tarball beside it.

    A `terminate` with its own `basename` and a different `package: basename:` writes the NAMD
    config under the package basename but ships the state files -- the minimal .prm among them --
    under the terminate basename.  Re-deriving the consolidated parameter file's name from the
    basename in force made the config reference `{package}_minimal.prm` while the tarball carried
    `{terminate}_minimal.prm`, so unpacking and running the package -- the whole point of packaging
    -- died on a missing parameter file.  The mismatch is invisible in the build directory, where
    both names exist.
    """

    def test_every_file_named_in_the_config_is_packaged(self):
        config = Config().configure_new()
        base, pkg = 'mysys', 'prod_mysys'
        with tempfile.TemporaryDirectory() as tmp:
            cwd = os.getcwd()
            os.chdir(tmp)
            try:
                shutil.copy(PSF, f'{base}.psf')
                for ext in ('pdb', 'coor', 'vel'):
                    open(f'{base}.{ext}', 'w').write('dummy\n')
                open(f'{base}_minimal.prm', 'w').write('* minimal\n*\n\nEND\n')

                task = TerminateTask(
                    specs=dict(basename=base,
                               package=dict(basename=pkg,
                                            namd=dict(ensemble='minimize', temperature=310,
                                                      nsteps=0, dcdfreq=0, xstfreq=0,
                                                      minimize=100))),
                    taskname='terminate', index=0)
                task.provision(dict(controller_index=0, pipeline=PipelineContext(),
                                    resource_manager=config.RM, scripters=config.scripters,
                                    namd_global_config=config['user']['namd']))
                state = task.register(dict(psf=PSFFileArtifact(f'{base}.psf'),
                                           pdb=PDBFileArtifact(f'{base}.pdb'),
                                           coor=NAMDCoorFileArtifact(f'{base}.coor'),
                                           xsc=None,
                                           vel=NAMDVelFileArtifact(f'{base}.vel')),
                                      key='state', artifact_type=StateArtifacts)
                state.minimal_prm = CharmmffParFileArtifact(data=f'{base}_minimal.prm', keep=True)
                state.data['minimal_prm'] = state.minimal_prm

                task.make_package()

                with tarfile.open(f'{pkg}.tar.gz') as tf:
                    packaged = {os.path.basename(n) for n in tf.getnames()}
                    cfgname = next(n for n in tf.getnames() if n.endswith('.namd'))
                    cfg = tf.extractfile(cfgname).read().decode()
            finally:
                os.chdir(cwd)

        self.assertEqual(os.path.basename(cfgname), f'{pkg}.namd')
        self.assertIn(f'{base}_minimal.prm', packaged)
        referenced = {}
        for line in cfg.splitlines():
            fields = line.split()
            if len(fields) == 2 and fields[0].lower() in NAMD_FILE_KEYS:
                referenced[fields[0]] = fields[1]
        # the config must actually name the files, not just fail to name missing ones
        self.assertIn('parameters', referenced)
        self.assertEqual(referenced['parameters'], f'{base}_minimal.prm')
        missing = {k: v for k, v in referenced.items() if os.path.basename(v) not in packaged}
        self.assertEqual(missing, {},
                         f'packaged {cfgname} references file(s) absent from the tarball: {missing}')


if __name__ == '__main__':
    unittest.main()


#: A dihedral quartet the shipped feb26 default set defines TWICE, with different force constants:
#: `par_all36_cgenff.prm` gives Kchi=1.25 and `toppar_all36_carb_imlab.str` appends 3.1 (derived by
#: analogy, under a `!DNAP` label).  `merge` is last-wins, so whichever file is read last decides.
CONFLICTED = 'NG2O1  CG2R61 CG2R61 NG2S3'

_PRM = """* synthetic parameter file standing in for par_all36_cgenff.prm
*

ATOMS
MASS  -1  NG2O1   14.00700
MASS  -1  CG2R61  12.01100
MASS  -1  NG2S3   14.00700

DIHEDRALS
{quartet}     1.2500   2   180.00 ! the .prm value

NONBONDED
NG2O1    0.0  -0.2000  1.8500
CG2R61   0.0  -0.0700  1.9924
NG2S3    0.0  -0.2000  1.8500

END
""".format(quartet=CONFLICTED)

_STR = """* synthetic stream file standing in for toppar_all36_carb_imlab.str
*

read param card flex append
* appended parameters
*

DIHEDRALS
{quartet}     3.1000   2   180.00 ! the .str value, appended later

END
""".format(quartet=CONFLICTED)


class TestPackagedParametersMatchWhatWasSimulated(unittest.TestCase):
    """The consolidated ``.prm`` must resolve a duplicated term the way the NAMD run resolved it.

    ``NAMDScripter`` loads ``standard['prm'] + standard['str']``; with last-wins merging the
    stream's value is the one the system was simulated with.  ``TerminateTask`` used to merge its
    registered artifacts first and append the standard set afterwards, so on a build with no MD
    step -- where the standard files were never staged as artifacts -- every ``.prm`` landed after
    the streams and the packaged file disagreed with the simulation.

    No example in the suite reaches that path (configs with a ``terminate`` and no MD task: 0),
    which is why this survived; hence a synthetic one here.
    """

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.prm = os.path.join(self.d, 'standard.prm')
        self.strf = os.path.join(self.d, 'appended.str')
        with open(self.prm, 'w') as fh:
            fh.write(_PRM)
        with open(self.strf, 'w') as fh:
            fh.write(_STR)
        self.cwd = os.getcwd()
        os.chdir(self.d)

    def tearDown(self):
        os.chdir(self.cwd)
        shutil.rmtree(self.d, ignore_errors=True)

    def _run_md_less_terminate(self):
        """A build that ran no MD: the stream is a registered artifact, the .prm is not."""
        task = mock.Mock(spec=TerminateTask)
        task.basename = 'pkg'
        task.build_stamp.return_value = 'test'
        state = mock.Mock()
        state.psf.name = PSF
        state.psf.exists.return_value = True

        def _artifact(which):
            if which == 'state':
                return state
            if which == 'charmmff_streamfiles':
                return [_FA(self.strf)]
            return None        # no parfile artifacts: nothing staged them

        task.get_current_artifact.side_effect = _artifact
        scripter = mock.Mock()
        scripter.fetch_standard_charmm_parameters.return_value = [self.prm]
        task.get_scripter.return_value = scripter

        # The consolidated file keeps only records whose atom types are in the PSF, so the PSF has
        # to contain the quartet or the comparison has nothing to compare.
        psf = mock.Mock()
        psf.atoms = [mock.Mock(atomtype=t) for t in CONFLICTED.split()]
        with mock.patch('pestifer.tasks.terminate.PSFContents', return_value=psf):
            out = TerminateTask.generate_minimal_params(task)
        self.assertIsNotNone(out, 'POSITIVE CONTROL: nothing was written, so nothing was checked')
        return open(out).read()

    def _kchi(self, text):
        for line in text.splitlines():
            f = line.split('!')[0].split()
            if len(f) >= 5 and f[:4] == CONFLICTED.split():
                return float(f[4])
        return None

    def test_the_stream_value_survives_as_it_does_in_a_namd_run(self):
        kchi = self._kchi(self._run_md_less_terminate())
        self.assertIsNotNone(kchi, 'POSITIVE CONTROL: the quartet is absent, so order proves nothing')
        self.assertEqual(kchi, 3.1,
                         'the packaged file must carry the value a NAMD run would use (streams '
                         'merged last), not the .prm value the old artifact-first order gave')

    def test_the_two_synthetic_files_really_do_disagree(self):
        """Guards the fixture: if both files carried the same constant, the test above could not
        fail no matter which order won."""
        self.assertEqual(self._kchi(_PRM), 1.25)
        self.assertEqual(self._kchi(_STR), 3.1)
        self.assertNotEqual(self._kchi(_PRM), self._kchi(_STR))
