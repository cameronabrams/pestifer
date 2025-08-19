
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