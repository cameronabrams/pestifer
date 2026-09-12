# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Small-molecule model compounds are ligands, whichever force-field stream defines them.

pestifer derives a residue's segtype from the FILE that defines it.  Acetate (ACET) is defined in
``toppar_all36_prot_model.str`` and acetone (ACO) in ``toppar_all36_carb_model.str``, so they
derived as ``protein`` and ``glycan``.  Example 5 carries acetate as a crystallographic ligand and
example 24 solvates in acetone -- 1,753 acetone molecules that were classified as sugar.  Both are
now curated as ``ligand``, matching the other organic solvents (DMSO, ACN).
"""
import unittest

from pestifer.core.labels import Labels, _load_derived_segtypes


SMALL_MOLECULES = {'ACET': 'acetate', 'ACO': 'acetone'}


class TestSmallMoleculeModelCompoundsAreLigands(unittest.TestCase):

    def test_each_is_classified_ligand(self):
        for r, desc in SMALL_MOLECULES.items():
            with self.subTest(residue=r):
                self.assertEqual(Labels.segtype_of_resname.get(r), 'ligand',
                                 f'{r} ({desc}) must be a ligand')

    def test_they_are_curated_not_derived(self):
        curated = Labels.segtypes['ligand']['resnames']
        for r in SMALL_MOLECULES:
            self.assertIn(r, curated, f'{r} must be curated, or the derivation reclaims it')

    def test_the_derivation_alone_would_still_get_them_wrong(self):
        # Guards the test above against passing vacuously: if the generated table ever stops
        # misclassifying these, the curated entries are no longer what makes them ligands.
        derived = _load_derived_segtypes()
        where = {r: st for st, names in derived.items() for r in names if r in SMALL_MOLECULES}
        self.assertEqual(where.get('ACET'), 'protein')
        self.assertEqual(where.get('ACO'), 'glycan')

    def test_consistent_with_the_other_organic_solvents(self):
        for r in ('DMSO', 'ACN'):
            self.assertEqual(Labels.segtype_of_resname.get(r), 'ligand')


if __name__ == '__main__':
    unittest.main()
