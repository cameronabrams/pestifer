import os
import shutil
import tarfile
import tempfile
import unittest

from pestifer.core.artifacts import (CharmmffParFileArtifact, NAMDCoorFileArtifact,
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
