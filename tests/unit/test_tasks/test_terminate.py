import os
import unittest

from pestifer.tasks.terminate import TerminateTask

PSF = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'inputs', 'existing.psf'))


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


if __name__ == '__main__':
    unittest.main()
