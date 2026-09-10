# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
``patches:`` hands a raw name to psfgen, which ignores one it does not recognize.  A typo
therefore costs a modification silently -- the build succeeds and the system simply lacks the
change that was asked for.  The force field is already indexed; these pin that it is consulted.
"""
import unittest

from pestifer.core.errors import PestiferBuildError
from pestifer.core.resourcemanager import ResourceManager
from pestifer.tasks.psfgen import PsfgenTask


class TestPatchNameValidation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.RM = ResourceManager()

    def _task(self, patches):
        t = PsfgenTask.__new__(PsfgenTask)
        t.resource_manager = self.RM
        t.specs = {'mods': {'patches': patches}}
        return t

    def test_a_real_patch_is_accepted(self):
        self._task(['ASPP:C:199'])._validate_patch_names()

    def test_several_real_patches_are_accepted(self):
        self._task(['ASPP:C:199', 'KAC:A:5'])._validate_patch_names()

    def test_a_typo_is_rejected_with_near_misses(self):
        with self.assertRaises(PestiferBuildError) as cm:
            self._task(['ASPPP:C:199'])._validate_patch_names()
        self.assertIn('ASPP', str(cm.exception))

    def test_a_residue_name_is_rejected_and_explained(self):
        # SEP/TPO/PTR/TYS ship as whole RESI residues; reaching for them here is the likely
        # real mistake, and "not found" would be a misleading thing to say about a name the
        # force field does define
        with self.assertRaises(PestiferBuildError) as cm:
            self._task(['SEP:A:12'])._validate_patch_names()
        msg = str(cm.exception)
        self.assertIn('RESI', msg)
        self.assertIn('mutation', msg)

    def test_no_patches_is_not_an_error(self):
        self._task([])._validate_patch_names()
        t = PsfgenTask.__new__(PsfgenTask)
        t.resource_manager = self.RM
        t.specs = {}
        t._validate_patch_names()

    def test_an_unavailable_index_is_not_treated_as_a_typo(self):
        t = PsfgenTask.__new__(PsfgenTask)
        t.resource_manager = None            # no RM at all
        t.specs = {'mods': {'patches': ['NONSENSE:A:1']}}
        t._validate_patch_names()            # must not raise
