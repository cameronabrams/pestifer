import importlib.util
import shutil
import tempfile
import unittest
from pathlib import Path

from pestifer.core import labels


class TestRegisterResnameSegtype(unittest.TestCase):
    """Source-editing of ``_segtypes`` in labels.py used by ``modify-package``."""

    def setUp(self):
        self.d = Path(tempfile.mkdtemp())
        self.tmp = self.d / 'labels.py'
        shutil.copy(labels.__file__, self.tmp)

    def _load(self):
        spec = importlib.util.spec_from_file_location('editedlabels', self.tmp)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        return m

    def test_insert_registers_and_classifies(self):
        path, added = labels.register_resnames_segtype(['MYLIG', 'myl2'], 'ligand', labels_path=self.tmp)
        self.assertEqual(added, ['MYLIG', 'MYL2'])   # upper-cased
        m = self._load()
        self.assertIn('MYLIG', m._segtypes['ligand']['resnames'])
        self.assertIn('MYL2', m._segtypes['ligand']['resnames'])
        # the segtype map derived at import time picks it up
        self.assertEqual(m.Labels.segtype_of_resname['MYLIG'], 'ligand')
        # existing names are preserved
        self.assertIn('LF0', m._segtypes['ligand']['resnames'])

    def test_idempotent_for_existing_name(self):
        _, added = labels.register_resnames_segtype('LF0', 'ligand', labels_path=self.tmp)
        self.assertEqual(added, [])

    def test_unknown_segtype_raises(self):
        with self.assertRaises(ValueError):
            labels.register_resnames_segtype('MYLIG', 'nonesuch', labels_path=self.tmp)

    def test_result_still_parses(self):
        labels.register_resnames_segtype('NEWION', 'ion', labels_path=self.tmp)
        # loading the edited module is itself a parse+exec check
        m = self._load()
        self.assertIn('NEWION', m._segtypes['ion']['resnames'])

    def test_segtype_names_lists_keys(self):
        names = labels.segtype_names()
        self.assertIn('ligand', names)
        self.assertIn('ion', names)


if __name__ == '__main__':
    unittest.main()
