import unittest

from pestifer.sphinxext.task_table import _single_task_details as detail


class TestSolvateDetail(unittest.TestCase):
    def test_water_is_default(self):
        self.assertEqual(detail('solvate', {}), 'water box')
        self.assertEqual(detail('solvate', {'solvent': 'TIP3'}), 'water box')
        self.assertEqual(detail('solvate', {'solvent': 'water'}), 'water box')

    def test_non_water_uses_solvent_name(self):
        self.assertEqual(detail('solvate', {'solvent': 'DMSO'}), 'DMSO box')
        self.assertEqual(detail('solvate', {'solvent': 'MEOH'}), 'MEOH box')

    def test_salt_and_ions(self):
        self.assertEqual(detail('solvate', {'solvent': 'DMSO', 'salt_con': 0.15}),
                         'DMSO box + 0.15 M salt')
        self.assertEqual(detail('solvate', {'salt_con': 0.15, 'cation': 'SOD', 'anion': 'CLA'}),
                         'water box + 0.15 M SOD/CLA')


if __name__ == '__main__':
    unittest.main()
