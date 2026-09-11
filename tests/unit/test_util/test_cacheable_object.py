
import tarfile
import yaml
import os
from pestifer.util.cacheable_object import CacheableObject, TarBytesFS
import unittest
import random
from pathlib import Path
# define a subclass of CacheableObject
class CacheableObjectSubclass(CacheableObject):
    def _build_from_resources(self, resource_root: Path) -> None:
        # read in mock_db.yaml
        with open(os.path.join(resource_root, 'mock_db.yaml'), 'r') as f:
            data = yaml.safe_load(f)
        self.data = data
        self.from_cache = False


class TestTarBytesFS(unittest.TestCase):

    def setUp(self):
        if not os.path.isfile('resources/mock_db.yaml'):
            # make a mock database in dictionary form with random three-letter keys and random floats for values
            mock_db = {f'key{i}': random.random() for i in range(1, 6)}
            # write this as YAML to resources
            with open('resources/mock_db.yaml', 'w') as f:
                yaml.dump(mock_db, f)
            # now tar it
        if not os.path.isfile('resources/mock_db.tar.gz'):
            with tarfile.open('resources/mock_db.tar.gz', 'w:gz') as tar:
                tar.add('resources/mock_db.yaml', arcname='mock_db.yaml')
        # and save the tar to disk
            
    def test_ls(self):
        # Test listing files in the tar
        self.tar_fs = TarBytesFS.from_file('resources/mock_db.tar.gz', compression='gzip')
        files = self.tar_fs.ls()
        self.assertIn('mock_db.yaml', [x['name'] for x in files])

    def test_open(self):
        # Test opening a file in the tar
        self.tar_fs = TarBytesFS.from_file('resources/mock_db.tar.gz', compression='gzip')
        with self.tar_fs.open('mock_db.yaml', 'r') as f:
            data = f.read()
            self.assertIn('key1', data)

class TestCacheableObject(unittest.TestCase):

    def setUp(self):
        if not os.path.isfile('resources/mock_db.yaml'):
            # make a mock database in dictionary form with random three-letter keys and random floats for values
            mock_db = {f'key{i}': random.random() for i in range(1, 6)}
            # write this as YAML to resources
            with open('resources/mock_db.yaml', 'w') as f:
                yaml.dump(mock_db, f)

    def test_cacheable_object(self):

        obj = CacheableObjectSubclass('resources', force_rebuild=True)
        self.assertIsInstance(obj, CacheableObjectSubclass)
        self.assertFalse(obj.from_cache)

        another_obj = CacheableObjectSubclass('resources')
        self.assertIsInstance(another_obj, CacheableObjectSubclass)
        self.assertTrue(another_obj.from_cache)

        self.assertEqual(obj.data, another_obj.data)

class TestCacheManagement(unittest.TestCase):
    def test_cache_files_and_clear(self):
        import tempfile
        from unittest import mock
        with tempfile.TemporaryDirectory() as d:
            dpath = Path(d)
            (dpath / 'cacheobj-foo-abc-v2.8.joblib').write_bytes(b'x')
            (dpath / 'cacheobj-bar-def-v2.8.joblib').write_bytes(b'y')
            (dpath / 'cacheobj-foo-abc-v2.8.joblib.lock').write_bytes(b'')
            (dpath / 'unrelated.txt').write_text('keep me')
            with mock.patch.object(CacheableObject, 'cache_directory', classmethod(lambda cls: dpath)):
                files = CacheableObject.cache_files()
                self.assertEqual(len(files), 2)
                self.assertTrue(all(f.suffix == '.joblib' for f in files))
                removed = CacheableObject.clear_cache()
                self.assertEqual(len(removed), 2)
                self.assertEqual(CacheableObject.cache_files(), [])
            # non-cache files are left alone; lock files are removed
            self.assertTrue((dpath / 'unrelated.txt').exists())
            self.assertFalse((dpath / 'cacheobj-foo-abc-v2.8.joblib.lock').exists())


class TestCacheKeyIsPerResourceRoot(unittest.TestCase):
    """
    A caller-supplied ``resource_label`` names the force-field RELEASE, which two pestifer
    installations share while their files do not -- the cached object stores absolute paths,
    including to the package's own ``charmmff/custom/`` files beside ``resource_root``.  Keyed
    on the label alone, a run from a throwaway tree (a clean-export CI check, a tox env, a pip
    install in a container) overwrote the shared entry with paths into itself, and every later
    run from the real install died on a FileNotFoundError naming a directory that no longer
    existed.
    """

    def _root(self, tmp, name, payload):
        import yaml as _yaml
        d = Path(tmp) / name
        d.mkdir(parents=True)
        with open(d / 'mock_db.yaml', 'w') as f:
            _yaml.dump(payload, f)
        return d

    def test_same_label_different_roots_do_not_share_a_cache_entry(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            cdir = Path(tmp) / 'cache'
            a = self._root(tmp, 'install_a/feb26', {'who': 'A'})
            b = self._root(tmp, 'install_b/feb26', {'who': 'B'})

            oa = CacheableObjectSubclass(a, cache_dir=cdir, resource_label='feb26')
            ob = CacheableObjectSubclass(b, cache_dir=cdir, resource_label='feb26')
            self.assertEqual(oa.data['who'], 'A')
            # the direction that actually breaks: with one shared entry the SECOND install is
            # served the first's cache, because the entry is newer than its own resources
            self.assertEqual(ob.data['who'], 'B', "B was served A's cache")

            files = sorted(p.name for p in cdir.glob('*.joblib'))
            self.assertEqual(len(files), 2, f'expected one cache entry per root, got {files}')
            # the human-readable release name survives in both
            self.assertTrue(all('feb26' in f for f in files), files)

    def test_one_root_still_reuses_its_entry(self):
        # A guard against over-fixing, not a bug catcher: this passes with or without the key
        # change, and exists so that disambiguating the key cannot silently disable caching.
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            cdir = Path(tmp) / 'cache'
            a = self._root(tmp, 'install_a/feb26', {'who': 'A'})
            first = CacheableObjectSubclass(a, cache_dir=cdir, resource_label='feb26')
            second = CacheableObjectSubclass(a, cache_dir=cdir, resource_label='feb26')
            self.assertFalse(first.from_cache)
            self.assertTrue(second.from_cache, 'caching must still work for a single root')
