# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import os
import signal
import subprocess
import time
import unittest

from pestifer.core.command import Command
from pestifer.core import command as command_mod


class TestCommand(unittest.TestCase):
    def test_echo(self):
        expected_stdout='123\n'
        c=Command('echo 123')
        c.run()
        self.assertEqual(expected_stdout,c.stdout)

    def test_run_unregisters_child_when_done(self):
        # a completed command must leave nothing behind in the active-child registry
        before = set(command_mod._active_child_pgids)
        Command('echo hi').run(quiet=True)
        self.assertEqual(set(command_mod._active_child_pgids), before)


class TestChildShutdown(unittest.TestCase):
    """The teardown machinery that turns Ctrl-C / kill into a clean shutdown of the
    whole external-process tree (charmrun -> namd3 -> PEs)."""

    @staticmethod
    def _pgroup_members(pgid):
        out = subprocess.run(['pgrep', '-g', str(pgid)], capture_output=True, text=True).stdout
        return [int(x) for x in out.split()]

    def test_terminate_children_kills_whole_group(self):
        # a shell that backgrounds a sleep and waits -> a 2-process group (sh + sleep)
        p = subprocess.Popen('sleep 60 & wait', shell=True, start_new_session=True)
        pgid = os.getpgid(p.pid)
        self.addCleanup(lambda: self._pgroup_members(pgid) and os.killpg(pgid, signal.SIGKILL))
        time.sleep(0.4)
        self.assertGreaterEqual(len(self._pgroup_members(pgid)), 2)   # sh + sleep both alive

        command_mod._register_child(p.pid)
        n = command_mod._terminate_children(signal.SIGKILL)
        command_mod._unregister_child(p.pid)
        self.assertEqual(n, 1)                                        # one group signaled
        try:
            p.wait(timeout=3)
        except subprocess.TimeoutExpired:
            pass
        time.sleep(0.3)
        self.assertEqual(self._pgroup_members(pgid), [],              # whole group gone
                         'process group survived the teardown')

    def test_register_unregister(self):
        command_mod._register_child(999999)
        self.assertIn(999999, command_mod._active_child_pgids)
        command_mod._unregister_child(999999)
        self.assertNotIn(999999, command_mod._active_child_pgids)
        # terminate on a dead/absent pid must not raise
        command_mod._register_child(999999)
        command_mod._terminate_children(signal.SIGTERM)
        command_mod._unregister_child(999999)

    def test_install_signal_handlers_is_noop_under_pytest(self):
        # under pytest, install must leave the test runner's own SIGINT handling intact
        before = signal.getsignal(signal.SIGINT)
        command_mod.install_signal_handlers()
        self.assertIs(signal.getsignal(signal.SIGINT), before)
        self.assertIsNot(signal.getsignal(signal.SIGINT), command_mod._signal_handler)
