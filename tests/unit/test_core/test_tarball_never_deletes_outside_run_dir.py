# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A build archives its intermediates; it must never delete a file it did not create.

``FileArtifactList.make_tarball(remove=True)`` archives each artifact and then ``os.remove``s it.
There was no check that the path lay inside the run directory, and ``terminate`` calls it that
way twice (``terminate.py:232`` and ``:402``). Any artifact registered from a user-supplied path
was therefore deleted from wherever it lived.

Reported 2026-10-12: a ``bilayer.prebuilt`` bilayer carrying **28 ns of GPU calibration**, in a
sibling run's folder, was archived into a build's tarball and removed from its own directory. The
files it did not name (``.coor``/``.vel``) stayed behind, so the source folder still looked
populated — which is the shape that makes this expensive to notice.

Two fixes, because one of them alone leaves the next caller exposed:

* ``make_membrane_system``'s prebuilt branch copies its inputs in and registers the COPIES, which
  is what ``continuation`` has always done for its five inputs.
* ``make_tarball`` refuses to remove anything outside the run directory, whatever registered it.

The second is the one that generalizes. ``_ensure_in_cwd`` lived in ``continuation.py``, private,
called nowhere else — which is precisely why the second consumer of that idea was written without
it. It is shared now (:func:`pestifer.core.artifacts.ensure_in_cwd`), but a guard at the
destructive operation is what protects a consumer nobody has written yet.
"""
import os
import tarfile
import tempfile
import unittest
from pathlib import Path

from pestifer.core.artifacts import (FileArtifactList, TXTFileArtifact, _inside_cwd,
                                     ensure_in_cwd)


class TestInsideCwd(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.outside = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.tmp)

    def tearDown(self):
        os.chdir(self.cwd)

    def test_a_file_in_the_run_dir_is_inside(self):
        Path('a.txt').write_text('x')
        self.assertTrue(_inside_cwd('a.txt'))
        self.assertTrue(_inside_cwd(Path('a.txt').resolve()))

    def test_a_file_in_a_subdirectory_is_inside(self):
        Path('sub').mkdir()
        Path('sub/a.txt').write_text('x')
        self.assertTrue(_inside_cwd('sub/a.txt'))

    def test_a_sibling_directory_is_OUTSIDE(self):
        p = Path(self.outside) / 'precious.psf'
        p.write_text('28 ns of GPU time')
        self.assertFalse(_inside_cwd(p))

    def test_dotdot_cannot_walk_out(self):
        """NEGATIVE CONTROL for a string-prefix implementation."""
        p = Path(self.outside) / 'precious.psf'
        p.write_text('x')
        rel = os.path.relpath(p, start=self.tmp)       # '../<outside>/precious.psf'
        self.assertTrue(rel.startswith('..'), f'control is not a dotdot path: {rel!r}')
        self.assertFalse(_inside_cwd(rel))
        self.assertFalse(_inside_cwd(os.path.join('sub', '..', rel)))

    def test_a_symlink_IN_the_run_dir_counts_as_inside(self):
        """Because `os.remove` unlinks the entry, not what it points at.

        Treating such a link as outside is not "safer" -- it cannot protect the target, which was
        never at risk, and it leaves a stray link behind every run. Several test fixtures are
        exactly this shape and did leave strays until the guard was narrowed to the entry.
        """
        target = Path(self.outside) / 'precious.psf'
        target.write_text('x')
        link = Path('looks-local.psf')
        try:
            link.symlink_to(target)
        except (OSError, NotImplementedError):
            self.skipTest('symlinks unavailable')
        self.assertTrue(_inside_cwd(link))
        os.remove(link)
        self.assertTrue(target.exists(), 'removing a symlink must not touch its target')

    def test_a_symlinked_DIRECTORY_cannot_launder_an_outside_file(self):
        """NEGATIVE CONTROL for narrowing the guard: the parent is still resolved."""
        target_dir = Path(self.outside)
        (target_dir / 'precious.psf').write_text('x')
        try:
            Path('linkdir').symlink_to(target_dir, target_is_directory=True)
        except (OSError, NotImplementedError):
            self.skipTest('symlinks unavailable')
        self.assertFalse(_inside_cwd('linkdir/precious.psf'))


class TestMakeTarballRemoval(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.outside = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.tmp)
        self.local = Path('intermediate.txt')
        self.local.write_text('local build intermediate')
        self.foreign = Path(self.outside) / 'equil2-final.txt'
        self.foreign.write_text('28 ns of GPU calibration')

    def tearDown(self):
        os.chdir(self.cwd)

    def _list(self):
        return FileArtifactList([TXTFileArtifact(str(self.local)),
                                 TXTFileArtifact(str(self.foreign))])

    def test_the_foreign_file_SURVIVES(self):
        """The bug, stated directly."""
        self._list().make_tarball('archive', remove=True)
        self.assertTrue(self.foreign.exists(),
                        'make_tarball deleted a file outside the run directory')

    def test_the_local_intermediate_is_still_removed(self):
        """POSITIVE CONTROL: the guard must not disable the feature it is guarding."""
        self._list().make_tarball('archive', remove=True)
        self.assertFalse(self.local.exists(),
                         'the guard stopped make_tarball removing its own build intermediates')

    def test_BOTH_are_archived(self):
        """Keeping the foreign file must not mean dropping it from the tarball."""
        self._list().make_tarball('archive', remove=True)
        with tarfile.open('archive.tar.gz') as tf:
            names = {Path(n).name for n in tf.getnames()}
        self.assertIn('intermediate.txt', names)
        self.assertIn('equil2-final.txt', names)

    def test_it_says_so(self):
        with self.assertLogs('pestifer.core.artifacts', level='WARNING') as cm:
            self._list().make_tarball('archive', remove=True)
        msg = '\n'.join(cm.output)
        self.assertIn('equil2-final.txt', msg)
        self.assertIn('NOT', msg, 'the warning does not say the file was kept')


class TestEnsureInCwd(unittest.TestCase):
    """The site fix: register the copy, not the caller's original."""

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.outside = tempfile.mkdtemp()
        self.cwd = os.getcwd()
        os.chdir(self.tmp)

    def tearDown(self):
        os.chdir(self.cwd)

    def test_an_outside_file_is_copied_in_and_the_basename_returned(self):
        src = Path(self.outside) / 'equil2-final.psf'
        src.write_text('precious')
        got = ensure_in_cwd(str(src))
        self.assertEqual(got, 'equil2-final.psf')
        self.assertTrue(Path(got).exists())
        self.assertTrue(src.exists(), 'the original must not be moved')
        self.assertFalse(_inside_cwd(src))
        self.assertTrue(_inside_cwd(got),
                        'the returned path must be the one it is safe to register')

    def test_a_local_file_is_left_alone(self):
        p = Path('already-here.psf')
        p.write_text('x')
        self.assertEqual(ensure_in_cwd('already-here.psf'), 'already-here.psf')

    def test_a_source_that_does_not_exist_is_passed_through(self):
        """Not a weakening of the fix, and the reason is structural.

        `make_tarball` archives and removes strictly inside `if artifact.exists()`, so a path
        with no file behind it can never be deleted -- there is nothing to protect. Copying would
        raise here and replace the caller's specific error ("prebuilt xsc missing") with a
        FileNotFoundError from a helper. Pinned so a later "tighten this" does not reintroduce
        a copy that fails on a path the caller is about to diagnose better.
        """
        missing = str(Path(self.outside) / 'nope.psf')
        self.assertEqual(ensure_in_cwd(missing), missing)
        self.assertFalse(Path('nope.psf').exists())

    def test_a_nonexistent_artifact_is_never_removed_anyway(self):
        """The structural claim above, verified rather than asserted."""
        gone = Path(self.outside) / 'vanished.txt'
        FileArtifactList([TXTFileArtifact(str(gone))]).make_tarball('a', remove=True)
        self.assertFalse(gone.exists())   # it never existed; the point is nothing raised

    def test_empty_is_passed_through(self):
        for empty in ('', None):
            self.assertEqual(ensure_in_cwd(empty), empty)


class TestThePrebuiltBranchUsesIt(unittest.TestCase):
    """The idea had one consumer and a second was written without it. Pin the second."""

    def test_the_prebuilt_branch_registers_copies(self):
        import inspect
        from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask
        src = inspect.getsource(MakeMembraneSystemTask)
        i = src.index("key='quilt_state'")
        window = src[max(0, i - 600):i]
        self.assertIn('ensure_in_cwd', window,
                      "the prebuilt branch registers the user's own paths again, so terminate "
                      "will delete their bilayer from its own folder")
        for k in ('pdb', 'psf', 'xsc'):
            self.assertIn(f"ensure_in_cwd(self.bilayer_specs['prebuilt']['{k}'])", window,
                          f'{k} is still registered unprotected')


if __name__ == '__main__':
    unittest.main()
