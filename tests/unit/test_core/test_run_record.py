# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the machine-readable record of what a build did."""

import json
import os
import tempfile
import unittest
from unittest import mock
from types import SimpleNamespace

from pestifer.core import run_record as RR


class _Task:
    def __init__(self, index, taskname, outcome=None, system_facts=None):
        self.index, self.taskname = index, taskname
        self.outcome = outcome or {}
        if system_facts is not None:
            self.system_facts = system_facts


class TestProtocol(unittest.TestCase):

    def test_only_tasks_that_ran_something_appear(self):
        tasks = [_Task(0, 'fetch'), _Task(1, 'psfgen'),
                 _Task(2, 'md', {'ensemble': 'minimize', 'steps': 100})]
        protocol = RR.protocol_from_tasks(tasks)
        self.assertEqual([p['task'] for p in protocol], ['md'])

    def test_order_is_preserved(self):
        tasks = [_Task(i, f't{i}', {'steps': i}) for i in range(4)]
        self.assertEqual([p['index'] for p in RR.protocol_from_tasks(tasks)], [0, 1, 2, 3])

    def test_adaptive_outcome_carries_its_stop_reason(self):
        tasks = [_Task(5, 'density_equilibrate',
                       {'adaptive': True, 'converged': False, 'steps': 8000,
                        'stopped_because': 'CEILING: ...'})]
        entry = RR.protocol_from_tasks(tasks)[0]
        self.assertFalse(entry['converged'])
        self.assertIn('CEILING', entry['stopped_because'])

    def test_no_tasks_is_not_an_error(self):
        self.assertEqual(RR.protocol_from_tasks(None), [])


class TestSystemFacts(unittest.TestCase):

    def test_taken_from_the_task_that_captured_them(self):
        # terminate reads them before packaging archives the files; they cannot be re-read later
        tasks = [_Task(0, 'psfgen'), _Task(9, 'terminate', system_facts={'atoms': 42})]
        self.assertEqual(RR.system_facts_from_tasks(tasks), {'atoms': 42})

    def test_absent_when_nothing_captured(self):
        self.assertEqual(RR.system_facts_from_tasks([_Task(0, 'psfgen')]), {})

    def test_missing_state_is_not_an_error(self):
        self.assertEqual(RR.system_facts(None), {})


class TestVersionGuard(unittest.TestCase):

    def _write(self, td, payload):
        path = os.path.join(td, RR.RUN_RECORD_NAME)
        with open(path, 'w') as f:
            json.dump(payload, f)
        return path

    def test_rejects_a_record_from_a_future_pestifer(self):
        # better to refuse than to misread a shape we do not know
        with tempfile.TemporaryDirectory() as td:
            path = self._write(td, {'run_record_version': RR.RUN_RECORD_VERSION + 1})
            with self.assertRaises(ValueError):
                RR.load_run_record(path)

    def test_rejects_a_file_that_is_not_a_record(self):
        with tempfile.TemporaryDirectory() as td:
            path = self._write(td, {'something': 'else'})
            with self.assertRaises(ValueError):
                RR.load_run_record(path)

    def test_accepts_the_current_version(self):
        with tempfile.TemporaryDirectory() as td:
            path = self._write(td, {'run_record_version': RR.RUN_RECORD_VERSION, 'seed': 1})
            self.assertEqual(RR.load_run_record(path)['seed'], 1)

    def test_load_accepts_directories_and_preserves_replica_order(self):
        with tempfile.TemporaryDirectory() as td:
            dirs = []
            for i in (1, 2, 3):
                d = os.path.join(td, f'rep-{i}')
                os.makedirs(d)
                with open(os.path.join(d, RR.RUN_RECORD_NAME), 'w') as f:
                    json.dump({'run_record_version': RR.RUN_RECORD_VERSION, 'seed': i}, f)
                dirs.append(d)
            self.assertEqual([r['seed'] for r in RR.load_run_records(dirs)], [1, 2, 3])


class TestBuildRecord(unittest.TestCase):

    def test_carries_the_reports_it_is_given(self):
        cfg = {'user': {'title': 'x', 'namd': {'seed': 7}}}
        record = RR.build_run_record(cfg, [], environment={'charmmff': 'feb26'},
                                     citations={'entries': []})
        self.assertEqual(record['seed'], 7)
        self.assertEqual(record['environment']['charmmff'], 'feb26')
        self.assertEqual(record['run_record_version'], RR.RUN_RECORD_VERSION)

    def test_survives_an_unreadable_config(self):
        record = RR.build_run_record(None, [])
        self.assertEqual(record['run_record_version'], RR.RUN_RECORD_VERSION)

    def test_write_never_raises(self):
        self.assertIsNone(RR.write_run_record({'a': 1}, path='/nonexistent-dir-xyzzy/rec.json'))


class TestCompletionIsStated(unittest.TestCase):
    """A record is only written for a build that finished, so its existence already implies
    success.  But *absence* is ambiguous -- failed, crashed, or still running -- and a sweep asking
    "did all 81 replicas succeed?" should be able to answer from a file's contents rather than from
    a directory listing.
    """

    def _record(self):
        cfg = mock.Mock()
        cfg.__getitem__ = lambda _s, k: {} if k == 'user' else None
        cfg.userfile = 'x.yaml'
        cfg.namd_type = 'cpu'
        return RR.build_run_record(cfg, [])

    def test_the_record_says_it_completed(self):
        rec = self._record()
        self.assertEqual(rec['status'], 'completed')
        self.assertEqual(rec['exit_code'], 0)

    def test_the_field_is_a_real_zero_not_a_missing_key(self):
        """`.get('exit_code')` returning None for an absent key reads exactly like a recorded
        failure, which is how this was misdiagnosed in the first place."""
        rec = self._record()
        self.assertIn('exit_code', rec)
        self.assertIsNotNone(rec.get('exit_code'))

    def test_the_record_version_moved_with_the_shape(self):
        self.assertGreaterEqual(RR.RUN_RECORD_VERSION, 2)
