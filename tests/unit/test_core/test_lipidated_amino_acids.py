# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
The six lipidated amino acids are amino acids, and must be classified as such.

pestifer derives a residue's segtype from which force-field FILE defines it.  That is right for
almost everything and wrong for these: all six live in ``stream/lipid/toppar_all36_lipid_prot.str``,
so they were derived as ``lipid`` and a palmitoylated cysteine could not be built into the
protein segment it belongs to.
"""
import unittest

from pestifer.core.labels import Labels


LIPIDATED = {
    'CYSP': 'S-palmitoyl-cysteine',
    'CYSF': 'S-farnesyl-cysteine',
    'CYSG': 'S-geranylgeranyl-cysteine',
    'CYSL': 'S-triacylhexadecane-cysteine',
    'GLYM': 'N-myristoyl-glycine',
    'LYSM': 'N-myristoyl-lysine',
}


class TestLipidatedAminoAcidsAreProtein(unittest.TestCase):

    def test_each_is_classified_protein(self):
        for r, desc in LIPIDATED.items():
            with self.subTest(residue=r):
                self.assertEqual(Labels.segtype_of_resname.get(r), 'protein',
                                 f'{r} ({desc}) must be protein, not lipid')

    def test_they_are_curated_not_derived(self):
        # the derivation keys on the defining file and gets these wrong, so they have to be in
        # the curated list -- which is what makes them win
        curated = Labels.segtypes['protein']['resnames']
        for r in LIPIDATED:
            self.assertIn(r, curated, f'{r} must be curated, or the derivation reclaims it')

    def test_real_lipids_are_still_lipids(self):
        # the fix must not reclassify the bilayer
        for r in ('POPC', 'POPE', 'CHL1', 'PSM', 'DPPC'):
            with self.subTest(residue=r):
                self.assertEqual(Labels.segtype_of_resname.get(r), 'lipid')

    def test_the_unmodified_parents_are_unaffected(self):
        for r in ('CYS', 'GLY', 'LYS'):
            self.assertEqual(Labels.segtype_of_resname.get(r), 'protein')
