import unittest

from pestifer.charmmff.segtype_classifier import segtype_of_topfile, derive_segtypes


class TestSegtypeClassifier(unittest.TestCase):
    def test_topfile_rules(self):
        # the systematic CHARMM file families map to segtypes
        self.assertEqual(segtype_of_topfile('top_all36_prot.rtf'), 'protein')
        self.assertEqual(segtype_of_topfile('top_all36_lipid.rtf'), 'lipid')
        self.assertEqual(segtype_of_topfile('toppar_all36_lipid_cholesterol.str'), 'lipid')
        self.assertEqual(segtype_of_topfile('top_all36_carb.rtf'), 'glycan')
        self.assertEqual(segtype_of_topfile('top_all36_na.rtf'), 'nucleicacid')
        self.assertEqual(segtype_of_topfile('top_all36_cgenff.rtf'), 'ligand')
        self.assertEqual(segtype_of_topfile('83G-cgenff.str'), 'ligand')
        self.assertEqual(segtype_of_topfile('toppar_all36_moreions.str'), 'ion')

    def test_none_when_unmatched(self):
        # a PDB-style alias not defined in a force-field file has no topfile
        self.assertIsNone(segtype_of_topfile(None))
        self.assertIsNone(segtype_of_topfile(''))
        self.assertIsNone(segtype_of_topfile('something_unrecognized.str'))

    def test_water_ions_split(self):
        # both water and ions live in water_ions.str; the water set decides
        water = {'TIP3', 'HOH'}
        self.assertEqual(segtype_of_topfile('toppar_water_ions.str', water_resnames=water, resname='TIP3'), 'water')
        self.assertEqual(segtype_of_topfile('toppar_water_ions.str', water_resnames=water, resname='SOD'), 'ion')

    def test_derive_excludes_curated_and_sorts(self):
        m = {
            'POPC': 'top_all36_lipid.rtf',
            'CHL1': 'toppar_all36_lipid_cholesterol.str',
            'AGLC': 'top_all36_carb.rtf',
            'ALA': 'top_all36_prot.rtf',
            'SOD': 'toppar_water_ions.str',
            'TIP3': 'toppar_water_ions.str',
            'NAG': None,   # PDB alias -> not derivable
        }
        derived = derive_segtypes(m, curated_names={'ALA'}, water_resnames={'TIP3'})
        self.assertEqual(derived['lipid'], ['CHL1', 'POPC'])   # sorted
        self.assertEqual(derived['glycan'], ['AGLC'])
        self.assertEqual(derived['ion'], ['SOD'])
        self.assertEqual(derived['water'], ['TIP3'])
        self.assertNotIn('protein', derived)   # ALA was curated -> excluded
        self.assertNotIn('NAG', [n for v in derived.values() for n in v])  # None topfile


if __name__ == '__main__':
    unittest.main()
