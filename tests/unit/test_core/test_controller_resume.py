# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""What a resumed run tells the rest of the build about itself.

Skipped tasks execute in no process and register nothing in this one's pipeline, so anything that
sweeps the directory or summarizes the build has to be told that the run resumed.  A 1.58M-atom
build's restart archived 26 files and left ~172 GB loose because nothing was.  These tests need no
toolchain: they drive ``do_tasks`` over fake tasks with the manifest stubbed out.
"""
import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest import mock

from pestifer.core.controller import Controller
from pestifer.core.run_manifest import MANIFEST_NAME, RunManifest


def _fake_task(index, taskname):
    t = SimpleNamespace(index=index, taskname=taskname, specs={}, result=0, duration=1.0,
                        run_resumed_from=None, executed=False)

    def execute():
        t.executed = True
        return 0

    t.execute = execute
    return t


class TestResumeIsVisibleToTasks(unittest.TestCase):

    def _run(self, resume_from):
        C = Controller.__new__(Controller)
        C.index = 0
        C.tasks = [_fake_task(i, n) for i, n in enumerate(['psfgen', 'md', 'terminate'])]
        C.pipeline = SimpleNamespace(get_current_artifact=lambda key: None)
        C.manifest, C.resume_from = None, 0
        manifest = mock.Mock()
        manifest.data = {'pestifer_version': '3.22.1', 'tasks': []}
        with mock.patch.object(Controller, '_init_run_manifest',
                               return_value=(manifest, resume_from)), \
             mock.patch.object(Controller, '_restore_state_for_resume', create=True), \
             mock.patch.object(Controller, '_clean_resumed_task_outputs', create=True):
            C.do_tasks()
        return C

    def test_every_task_learns_where_the_run_resumed(self):
        C = self._run(2)
        self.assertEqual([t.run_resumed_from for t in C.tasks], [2, 2, 2])
        self.assertEqual([t.executed for t in C.tasks], [False, False, True])

    def test_a_run_from_scratch_reports_no_resume(self):
        C = self._run(0)
        self.assertEqual([t.run_resumed_from for t in C.tasks], [0, 0, 0])
        self.assertTrue(all(t.executed for t in C.tasks))

    def test_the_manifest_and_resume_point_are_kept_for_the_run_record(self):
        # build.py assembles run-record.json from these after do_tasks returns
        C = self._run(2)
        self.assertEqual(C.resume_from, 2)
        self.assertEqual(C.manifest.data['pestifer_version'], '3.22.1')


class TestResumingKeepsTheBuildingVersion(unittest.TestCase):
    """``pestifer_version`` in the manifest -- and so in run-record.json -- must stay the version
    that built the system.  Overwriting it on resume let a six-minute terminate re-run claim a
    59-hour membrane build."""

    def _manifest_task(self, index, taskname, state):
        contract = SimpleNamespace(provides=frozenset(('state',)), requires=frozenset())
        return SimpleNamespace(index=index, taskname=taskname, specs={}, outcome={},
                               substage_outcomes=[], pipeline_contract=lambda s: contract)

    def test_an_existing_manifest_keeps_its_version_and_notes_the_resumer(self):
        with tempfile.TemporaryDirectory() as td:
            cwd = os.getcwd()
            os.chdir(td)
            try:
                psf = os.path.join(td, 'a.psf')
                open(psf, 'w').close()
                state = SimpleNamespace(psf=SimpleNamespace(name=psf), pdb=None, coor=None,
                                        xsc=None, vel=None)
                pipeline = SimpleNamespace(
                    get_current_artifact=lambda key: state if key == 'state' else None)
                built = RunManifest(os.path.join(td, MANIFEST_NAME), version='3.22.1')
                t0 = self._manifest_task(0, 'psfgen', state)
                built.record(t0, pipeline)

                C = Controller.__new__(Controller)
                C.index, C.restart, C.fresh, C.from_task = 0, True, False, None
                C.tasks = [t0, self._manifest_task(1, 'terminate', state)]
                manifest, resume_from = C._init_run_manifest()
            finally:
                os.chdir(cwd)

        self.assertEqual(resume_from, 1)                       # task 0 was already done
        from pestifer.util.stringthings import __pestifer_version__
        self.assertEqual(manifest.data['pestifer_version'], '3.22.1')
        self.assertEqual(manifest.data['resumed_by'], [__pestifer_version__])

