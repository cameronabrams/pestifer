# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the per-task run manifest (in-stream restart substrate)."""
import json
import os
import tempfile
import unittest
from types import SimpleNamespace

from pestifer.core.run_manifest import (
    RunManifest, spec_hash, tasks_fingerprint, MANIFEST_NAME,
)


def _fake_task(index, taskname, specs, provides=('STATE',)):
    contract = SimpleNamespace(provides=frozenset(provides))
    return SimpleNamespace(index=index, taskname=taskname, specs=specs,
                           pipeline_contract=lambda s: contract)


def _fake_pipeline(state_files):
    """A pipeline whose current 'state' exposes .psf/.pdb/... each with a .name."""
    attrs = {k: None for k in ('psf', 'pdb', 'coor', 'xsc', 'vel')}
    for k, v in state_files.items():
        attrs[k] = SimpleNamespace(name=v)
    state = SimpleNamespace(**attrs)
    return SimpleNamespace(get_current_artifact=lambda key: state if key == 'state' else None)


class TestSpecHash(unittest.TestCase):

    def test_order_independent(self):
        self.assertEqual(spec_hash({'a': 1, 'b': 2}), spec_hash({'b': 2, 'a': 1}))

    def test_sensitive_to_values(self):
        self.assertNotEqual(spec_hash({'nsteps': 100}), spec_hash({'nsteps': 200}))

    def test_handles_non_json_via_str(self):
        # default=str keeps it from raising on odd values
        spec_hash({'x': object()})   # should not raise


class TestTasksFingerprint(unittest.TestCase):

    def test_changes_on_spec_edit(self):
        a = [_fake_task(0, 'md', {'nsteps': 100})]
        b = [_fake_task(0, 'md', {'nsteps': 200})]
        self.assertNotEqual(tasks_fingerprint(a), tasks_fingerprint(b))

    def test_changes_on_task_insertion(self):
        a = [_fake_task(0, 'psfgen', {}), _fake_task(1, 'md', {})]
        b = [_fake_task(0, 'psfgen', {}), _fake_task(1, 'solvate', {}), _fake_task(2, 'md', {})]
        self.assertNotEqual(tasks_fingerprint(a), tasks_fingerprint(b))

    def test_stable_for_same_list(self):
        a = [_fake_task(0, 'psfgen', {'k': 1}), _fake_task(1, 'md', {'nsteps': 5})]
        b = [_fake_task(0, 'psfgen', {'k': 1}), _fake_task(1, 'md', {'nsteps': 5})]
        self.assertEqual(tasks_fingerprint(a), tasks_fingerprint(b))


class TestRunManifest(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.path = os.path.join(self.tmp, MANIFEST_NAME)

    def test_record_and_roundtrip(self):
        m = RunManifest(self.path, version='9.9.9', fingerprint='fp')
        m.record(_fake_task(0, 'psfgen', {}),
                 _fake_pipeline({'psf': '00-00-00_psfgen.psf', 'pdb': '00-00-00_psfgen.pdb'}))
        m.record(_fake_task(1, 'md', {'nsteps': 100}, provides=('STATE', 'MD_OUTPUT')),
                 _fake_pipeline({'psf': '00-00-00_psfgen.psf', 'pdb': '00-01-00_md.pdb',
                                 'coor': '00-01-00_md.coor', 'xsc': '00-01-00_md.xsc'}))
        # persisted and reloadable
        self.assertTrue(os.path.exists(self.path))
        r = RunManifest.load(self.path)
        self.assertEqual(r.data['pestifer_version'], '9.9.9')
        self.assertEqual(r.data['tasks_fingerprint'], 'fp')
        self.assertEqual([e['index'] for e in r.data['tasks']], [0, 1])
        md = r.data['tasks'][1]
        self.assertEqual(md['taskname'], 'md')
        self.assertEqual(md['state']['coor'], '00-01-00_md.coor')
        self.assertNotIn('vel', md['state'])   # absent slots omitted
        self.assertEqual(md['provides'], ['MD_OUTPUT', 'STATE'])

    def test_record_replaces_same_index(self):
        m = RunManifest(self.path)
        m.record(_fake_task(0, 'md', {'nsteps': 100}), _fake_pipeline({'psf': 'a.psf'}))
        m.record(_fake_task(0, 'md', {'nsteps': 200}), _fake_pipeline({'psf': 'b.psf'}))
        self.assertEqual(len(m.data['tasks']), 1)
        self.assertEqual(m.data['tasks'][0]['state']['psf'], 'b.psf')

    def test_mark_complete(self):
        m = RunManifest(self.path)
        m.record(_fake_task(0, 'md', {}), _fake_pipeline({'psf': 'a.psf'}))
        self.assertFalse(m.data['complete'])
        m.mark_complete()
        self.assertTrue(RunManifest.load(self.path).data['complete'])

    def test_record_never_raises(self):
        # a broken pipeline must not crash the build
        m = RunManifest(self.path)
        bad_pipeline = SimpleNamespace(get_current_artifact=lambda k: (_ for _ in ()).throw(RuntimeError('boom')))
        m.record(_fake_task(0, 'md', {}), bad_pipeline)   # swallowed + warned, no raise

    def test_save_is_atomic_no_partial(self):
        m = RunManifest(self.path)
        m.record(_fake_task(0, 'md', {}), _fake_pipeline({'psf': 'a.psf'}))
        # no leftover temp files in the dir
        leftovers = [f for f in os.listdir(self.tmp) if f.startswith('.manifest-')]
        self.assertEqual(leftovers, [])
        # and the file is valid JSON
        with open(self.path) as f:
            json.load(f)


class TestResumePoint(unittest.TestCase):
    """resume_point over a manifest x task-list x on-disk-file matrix."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.tmp)
        # a manifest for psfgen@0 -> md@1 -> solvate@2, with their state files touched on disk
        self.m = RunManifest(MANIFEST_NAME)
        self._record(0, 'psfgen', {}, {'psf': '00-00-00_psfgen.psf', 'pdb': '00-00-00_psfgen.pdb'})
        self._record(1, 'md', {'nsteps': 100}, {'psf': '00-00-00_psfgen.psf', 'pdb': '00-01-00_md.pdb',
                                                'coor': '00-01-00_md.coor'})
        self._record(2, 'solvate', {}, {'psf': '00-02-00_solvate.psf', 'pdb': '00-02-00_solvate.pdb'})

    def tearDown(self):
        os.chdir(self.cwd)

    def _record(self, i, name, specs, state):
        for fn in state.values():
            open(fn, 'w').close()   # touch
        self.m.record(_fake_task(i, name, specs), _fake_pipeline(state))

    def _tasks(self, *triples):
        return [_fake_task(i, n, s) for i, n, s in triples]

    def test_clean_full_match(self):
        tasks = self._tasks((0, 'psfgen', {}), (1, 'md', {'nsteps': 100}), (2, 'solvate', {}))
        self.assertEqual(self.m.resume_point(tasks), 2)   # all done -> resume at 3 (nothing)

    def test_spec_edit_invalidates_from_there(self):
        tasks = self._tasks((0, 'psfgen', {}), (1, 'md', {'nsteps': 999}), (2, 'solvate', {}))
        self.assertEqual(self.m.resume_point(tasks), 0)   # md changed -> resume at 1

    def test_taskname_change_breaks(self):
        tasks = self._tasks((0, 'psfgen', {}), (1, 'mdplot', {'nsteps': 100}))
        self.assertEqual(self.m.resume_point(tasks), 0)

    def test_missing_state_file_breaks(self):
        os.remove('00-01-00_md.coor')   # md's output vanished
        tasks = self._tasks((0, 'psfgen', {}), (1, 'md', {'nsteps': 100}), (2, 'solvate', {}))
        self.assertEqual(self.m.resume_point(tasks), 0)   # md no longer current -> resume at 1

    def test_complete_build_returns_minus_one(self):
        self.m.mark_complete()
        tasks = self._tasks((0, 'psfgen', {}), (1, 'md', {'nsteps': 100}), (2, 'solvate', {}))
        self.assertEqual(self.m.resume_point(tasks), -1)

    def test_no_manifest_entries_is_fresh(self):
        empty = RunManifest(MANIFEST_NAME + '.empty')
        self.assertEqual(empty.resume_point(self._tasks((0, 'psfgen', {}))), -1)

    def test_state_entry(self):
        self.assertEqual(self.m.state_entry(1)['coor'], '00-01-00_md.coor')
        self.assertEqual(self.m.state_entry(99), {})
