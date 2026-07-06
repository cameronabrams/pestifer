import os
import subprocess
import tempfile
import unittest
from pathlib import Path

from pestifer.util import gitutil


class TestGitContribute(unittest.TestCase):
    """The git helpers backing ``modify-package --branch``."""

    def setUp(self):
        self.d = Path(tempfile.mkdtemp())

        def g(*args):
            subprocess.run(['git', '-C', str(self.d), *args], check=True, capture_output=True)

        g('init', '-q', '-b', 'main')
        g('config', 'user.email', 't@t')
        g('config', 'user.name', 't')
        (self.d / 'a.txt').write_text('hello\n')
        g('add', '-A')
        g('commit', '-qm', 'init')

    def test_clean_and_current_branch(self):
        self.assertTrue(gitutil.worktree_is_clean(self.d))
        self.assertEqual(gitutil.current_branch(self.d), 'main')
        (self.d / 'a.txt').write_text('changed\n')
        self.assertFalse(gitutil.worktree_is_clean(self.d))

    def test_branch_exists(self):
        self.assertFalse(gitutil.branch_exists(self.d, 'nope'))
        subprocess.run(['git', '-C', str(self.d), 'branch', 'feature'], check=True, capture_output=True)
        self.assertTrue(gitutil.branch_exists(self.d, 'feature'))

    def test_create_branch_carries_changes_and_isolates(self):
        # working-tree changes (one tracked-modification, one new untracked file)
        # must be carried onto the new branch and committed there, leaving main clean
        (self.d / 'a.txt').write_text('modified\n')
        newf = self.d / 'new.txt'
        newf.write_text('contribution\n')
        gitutil.create_and_checkout_branch(self.d, 'contrib')
        self.assertEqual(gitutil.current_branch(self.d), 'contrib')
        gitutil.stage_and_commit(self.d, [str(self.d / 'a.txt'), str(newf)], 'contrib: change')
        self.assertTrue(gitutil.worktree_is_clean(self.d))
        subprocess.run(['git', '-C', str(self.d), 'checkout', '-q', 'main'], check=True, capture_output=True)
        # main has neither change
        self.assertEqual((self.d / 'a.txt').read_text(), 'hello\n')
        self.assertFalse((self.d / 'new.txt').exists())

    def test_create_existing_branch_raises(self):
        subprocess.run(['git', '-C', str(self.d), 'branch', 'dup'], check=True, capture_output=True)
        with self.assertRaises(gitutil.GitError):
            gitutil.create_and_checkout_branch(self.d, 'dup')

    def test_stage_and_commit_nothing_raises(self):
        # committing an unchanged path is an error, not a silent no-op
        with self.assertRaises(gitutil.GitError):
            gitutil.stage_and_commit(self.d, [str(self.d / 'a.txt')], 'noop')


if __name__ == '__main__':
    unittest.main()
