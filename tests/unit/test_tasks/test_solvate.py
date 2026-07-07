import unittest
from unittest import mock

from pestifer.tasks.solvate import SolvateTask
from pestifer.core.errors import PestiferError


def _make_task(repo):
    """A bare SolvateTask with just enough wiring to exercise _solvent_box_args."""
    t = SolvateTask.__new__(SolvateTask)
    cc = mock.Mock()
    cc.pdbrepository = repo
    rm = mock.Mock()
    rm.charmmff_content = cc
    t.resource_manager = rm
    return t


def _task_with_specs(specs):
    t = SolvateTask.__new__(SolvateTask)
    t.specs = specs
    return t


class TestIonizationPlan(unittest.TestCase):
    def test_water_always_ionizes(self):
        t = _task_with_specs({})
        self.assertEqual(t._ionization_plan('TIP3'), (True, True, False))
        self.assertEqual(t._ionization_plan(None)[0], True)

    def test_non_water_skips_ionization_by_default(self):
        do_ionize, is_water, ions_requested = _task_with_specs({})._ionization_plan('DMSO')
        self.assertFalse(do_ionize)
        self.assertFalse(is_water)
        self.assertFalse(ions_requested)

    def test_non_water_ionizes_when_ions_requested(self):
        for specs in ({'salt_con': 0.15}, {'cation': 'SOD'}, {'anion': 'CLA'}):
            do_ionize, is_water, ions_requested = _task_with_specs(specs)._ionization_plan('MEOH')
            self.assertTrue(do_ionize)
            self.assertTrue(ions_requested)
            self.assertFalse(is_water)


class TestSolventBoxArgs(unittest.TestCase):
    def test_water_uses_builtin_box(self):
        # every water alias returns no custom-box args (repo is never consulted)
        t = _make_task(repo=None)
        for w in ('TIP3', 'tip3', 'WATER', 'water', 'WAT', None):
            self.assertEqual(t._solvent_box_args(w), '')

    def test_non_water_box_entry(self):
        entry = mock.Mock()
        entry.is_box.return_value = True
        entry.get_box_psf.return_value = 'MEOH-box.psf'
        entry.get_box_pdb.return_value = 'MEOH-box.pdb'
        entry.get_box_edge.return_value = 20.34
        entry.get_key_atom.return_value = 'O1'
        repo = mock.MagicMock()
        repo.__contains__.return_value = True
        repo.checkout.return_value = entry
        t = _make_task(repo)
        args = t._solvent_box_args('MEOH')
        self.assertIn('-spsf MEOH-box.psf', args)
        self.assertIn('-spdb MEOH-box.pdb', args)
        self.assertIn('-ws 20.34', args)
        self.assertIn('-ks "name O1"', args)

    def test_non_water_missing_raises(self):
        repo = mock.MagicMock()
        repo.__contains__.return_value = False
        t = _make_task(repo)
        with self.assertRaises(PestiferError):
            t._solvent_box_args('NOPE')

    def test_non_water_molecule_entry_rejected(self):
        # a single-molecule entry (not a box) cannot be tiled by solvate
        entry = mock.Mock()
        entry.is_box.return_value = False
        repo = mock.MagicMock()
        repo.__contains__.return_value = True
        repo.checkout.return_value = entry
        t = _make_task(repo)
        with self.assertRaises(PestiferError):
            t._solvent_box_args('SOD')

    def test_key_atom_omitted_when_absent(self):
        entry = mock.Mock()
        entry.is_box.return_value = True
        entry.get_box_psf.return_value = 'X-box.psf'
        entry.get_box_pdb.return_value = 'X-box.pdb'
        entry.get_box_edge.return_value = 18.0
        entry.get_key_atom.return_value = None
        repo = mock.MagicMock()
        repo.__contains__.return_value = True
        repo.checkout.return_value = entry
        t = _make_task(repo)
        self.assertNotIn('-ks', t._solvent_box_args('X'))


if __name__ == '__main__':
    unittest.main()
