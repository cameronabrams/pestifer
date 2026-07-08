import argparse
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from pestifer.subcommands.modify_package import ModifyPackageSubcommand
from pestifer.core import ledger
from pestifer.util import gitutil


def _init_repo():
    d = tempfile.mkdtemp()
    subprocess.run(['git', 'init', '-q', d], check=True)
    subprocess.run(['git', '-C', d, 'config', 'user.email', 't@t'], check=True)
    subprocess.run(['git', '-C', d, 'config', 'user.name', 'tester'], check=True)
    (Path(d) / 'seed.txt').write_text('seed\n')
    subprocess.run(['git', '-C', d, 'add', '.'], check=True)
    subprocess.run(['git', '-C', d, 'commit', '-qm', 'init'], check=True)
    return Path(d)


def _parse(argv):
    sc = ModifyPackageSubcommand()
    p = argparse.ArgumentParser()
    subs = p.add_subparsers()
    sc.add_subparser(subs)
    return sc, p.parse_args(argv)


def _current_branch(d):
    return subprocess.run(['git', '-C', str(d), 'rev-parse', '--abbrev-ref', 'HEAD'],
                          capture_output=True, text=True).stdout.strip()


class TestModifyPackageBranchFlow(unittest.TestCase):
    def _run(self, d, argv, op_side_effect):
        sc, args = _parse(argv)
        rm = mock.Mock()
        rm.resources_path = str(d)      # the ledger lands inside the repo
        # route whichever example method the verb calls to the file-writing side effect
        for m in ('append_example', 'update_example', 'delete_example',
                  'rename_example', 'set_example_author'):
            getattr(rm, m).side_effect = op_side_effect
        with mock.patch('pestifer.subcommands.modify_package.ResourceManager', return_value=rm), \
             mock.patch.object(gitutil, 'package_repo_root', return_value=d):
            sc.func(args)

    def test_default_opens_branch_and_commits(self):
        d = _init_repo()
        start = _current_branch(d)

        def op(*a, **k):
            sub = d / 'examples' / '19'
            sub.mkdir(parents=True, exist_ok=True)
            (sub / 'x.rst').write_text('hi\n')

        self._run(d, ['modify-package', 'example', 'rename', '19', 'newname'], op)
        # a new auto-named contribution branch was created and checked out
        cur = _current_branch(d)
        self.assertNotEqual(cur, start)
        self.assertTrue(cur.startswith('modpkg/example-rename'), cur)
        # the change was committed (tree clean) and the file is in the tip commit
        self.assertTrue(gitutil.worktree_is_clean(d))
        show = subprocess.run(['git', '-C', str(d), 'show', '--name-only', '--format=', 'HEAD'],
                              capture_output=True, text=True).stdout
        self.assertIn('examples/19/x.rst', show)

    def test_modified_tracked_file_not_mangled(self):
        # a modified tracked file shows as " M path" (leading space) in git status; the
        # touched-path parser must not drop the first character of the path
        d = _init_repo()
        self._run(d, ['modify-package', 'example', 'rename', '19', 'n'],
                  lambda *a, **k: (d / 'seed.txt').write_text('changed\n'))
        self.assertTrue(gitutil.worktree_is_clean(d))          # committed, no mangled pathspec
        self.assertTrue(_current_branch(d).startswith('modpkg/example-rename'))
        show = subprocess.run(['git', '-C', str(d), 'show', '--name-only', '--format=', 'HEAD'],
                              capture_output=True, text=True).stdout
        self.assertIn('seed.txt', show)

    def test_explicit_branch_name(self):
        d = _init_repo()
        self._run(d, ['modify-package', 'example', 'add', 'foo.yaml', '--branch', 'my-feature'],
                  lambda *a, **k: (d / 'newfile.txt').write_text('x\n'))
        self.assertEqual(_current_branch(d), 'my-feature')

    def test_no_branch_applies_in_place(self):
        d = _init_repo()
        start = _current_branch(d)
        self._run(d, ['modify-package', 'example', 'rename', '19', 'n', '--no-branch'],
                  lambda *a, **k: (d / 'y.txt').write_text('x\n'))
        # stayed on the original branch, change left uncommitted
        self.assertEqual(_current_branch(d), start)
        self.assertFalse(gitutil.worktree_is_clean(d))

    def test_dirty_tree_rejected_by_default(self):
        d = _init_repo()
        (d / 'seed.txt').write_text('dirty\n')  # make the tree dirty
        sc, args = _parse(['modify-package', 'example', 'delete', '5'])
        with mock.patch('pestifer.subcommands.modify_package.ResourceManager', return_value=mock.Mock()), \
             mock.patch.object(gitutil, 'package_repo_root', return_value=d):
            with self.assertRaises(RuntimeError):
                sc.func(args)


class TestModifyPackageLedger(unittest.TestCase):
    def _rm(self, d):
        rm = mock.Mock()
        rm.resources_path = str(d)      # ledger lands inside the repo so it is committed
        for m in ('append_example', 'update_example', 'delete_example',
                  'rename_example', 'set_example_author'):
            getattr(rm, m).side_effect = lambda *a, **k: (d / 'foo.txt').write_text('hi\n')
        return rm

    def _run(self, d, rm, argv):
        sc, args = _parse(argv)
        with mock.patch('pestifer.subcommands.modify_package.ResourceManager', return_value=rm), \
             mock.patch.object(gitutil, 'package_repo_root', return_value=d):
            sc.func(args)

    def test_record_then_revert_roundtrip(self):
        d = _init_repo()
        rm = self._rm(d)

        # --- a modification is recorded in the ledger and committed ---
        self._run(d, rm, ['modify-package', 'example', 'rename', '19', 'newname'])
        entries = ledger.read(str(d))
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0]['category'], 'example')
        self.assertEqual(entries[0]['verb'], 'rename')
        self.assertIn('foo.txt', entries[0]['files'])
        self.assertFalse(entries[0]['reverted'])
        self.assertTrue((d / 'foo.txt').exists())
        self.assertTrue(gitutil.worktree_is_clean(d))       # ledger + change committed

        # --- reverting undoes the change and curates the ledger ---
        self._run(d, rm, ['modify-package', 'ledger', 'revert', '1'])
        self.assertFalse((d / 'foo.txt').exists())          # change reversed
        e1, e2 = ledger.get(str(d), 1), ledger.get(str(d), 2)
        self.assertTrue(e1['reverted'])                     # original marked reverted
        self.assertEqual(e1['reverted_by'], 2)
        self.assertEqual((e2['category'], e2['verb'], e2['reverts']), ('ledger', 'revert', 1))
        self.assertTrue(gitutil.worktree_is_clean(d))
        self.assertTrue(_current_branch(d).startswith('modpkg/revert-1'))

    def test_revert_already_reverted_is_rejected(self):
        d = _init_repo()
        rm = self._rm(d)
        self._run(d, rm, ['modify-package', 'example', 'rename', '19', 'n'])
        self._run(d, rm, ['modify-package', 'ledger', 'revert', '1'])
        with self.assertRaises(RuntimeError):
            self._run(d, rm, ['modify-package', 'ledger', 'revert', '1'])

    def test_revert_unknown_id_is_rejected(self):
        d = _init_repo()
        with self.assertRaises(RuntimeError):
            self._run(d, self._rm(d), ['modify-package', 'ledger', 'revert', '99'])

    def test_no_branch_records_without_committing(self):
        d = _init_repo()
        rm = self._rm(d)
        self._run(d, rm, ['modify-package', 'example', 'rename', '19', 'n', '--no-branch'])
        # recorded in the ledger, but left uncommitted in the working tree
        self.assertEqual(len(ledger.read(str(d))), 1)
        self.assertFalse(gitutil.worktree_is_clean(d))
        self.assertEqual(_current_branch(d), 'main' if gitutil.branch_exists(d, 'main') else 'master')


if __name__ == '__main__':
    unittest.main()
