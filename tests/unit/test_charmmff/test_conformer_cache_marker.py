# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
A conformer cache entry is a hit only when its commit finished.  ``shutil.move`` is atomic only
within one filesystem; across filesystems it degrades to copy-then-delete, so an interrupted
copy can leave a partial directory that an existence-only guard reads as a complete cache entry
and silently supplies as a conformer set.
"""
import inspect
import os
import tempfile
import unittest

from pestifer.charmmff.make_pdb_collection import (
    COMPLETION_MARKER, _cache_entry_is_complete, _mark_cache_entry_complete)


class TestConformerCacheCompleteness(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def _entry(self, name, files):
        d = os.path.join(self.tmp, name)
        os.mkdir(d)
        for f in files:
            open(os.path.join(d, f), 'w').write('x')
        return d

    def test_marked_entry_is_complete(self):
        d = self._entry('OK', ['a.pdb', 'a.psf'])
        _mark_cache_entry_complete(d)
        self.assertTrue(os.path.exists(os.path.join(d, COMPLETION_MARKER)))
        self.assertTrue(_cache_entry_is_complete(d))

    def test_partial_entry_is_rejected(self):
        # a copy interrupted before the PSF landed
        self.assertFalse(_cache_entry_is_complete(self._entry('PARTIAL', ['a.pdb'])))

    def test_empty_entry_is_rejected(self):
        self.assertFalse(_cache_entry_is_complete(self._entry('EMPTY', [])))

    def test_missing_entry_is_rejected(self):
        self.assertFalse(_cache_entry_is_complete(os.path.join(self.tmp, 'NOPE')))

    def test_a_file_is_not_an_entry(self):
        # sibling .lock files live beside entries and must not be mistaken for them
        p = os.path.join(self.tmp, '.X.lock')
        open(p, 'w').write('')
        self.assertFalse(_cache_entry_is_complete(p))

    def test_legacy_entry_without_marker_is_accepted(self):
        # caches written before the marker existed are intact; requiring the marker would force
        # every user to regenerate a working cache
        self.assertTrue(_cache_entry_is_complete(self._entry('LEGACY', ['a.pdb', 'a.psf'])))


class TestCompletenessCheckIsActuallyWired(unittest.TestCase):
    """The helper being correct is not enough -- the cache-hit guard has to consult it.  These
    tests exercise the helper directly, so without this the guard could be reverted to a bare
    ``os.path.exists`` and every test above would still pass."""

    def test_the_guard_consults_completeness_not_existence(self):
        from pestifer.charmmff import make_pdb_collection as m
        src = inspect.getsource(m)
        guard = [l for l in src.splitlines()
                 if 'successdir' in l and l.strip().startswith('if (')]
        self.assertEqual(len(guard), 1, f'expected one cache-hit guard, found {guard}')
        self.assertIn('_cache_entry_is_complete', guard[0])
        self.assertNotIn('os.path.exists', guard[0])
