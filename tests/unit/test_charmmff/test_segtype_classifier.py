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

    def test_model_compound_streams_are_ligands(self):
        # The four *_model.str streams hold small standalone molecules used to parameterize the
        # family their file is named for; matched by family they became protein/glycan/lipid/NA.
        for f in ('toppar_all36_prot_model.str', 'toppar_all36_carb_model.str',
                  'toppar_all36_lipid_model.str', 'toppar_all36_na_model.str'):
            with self.subTest(topfile=f):
                self.assertEqual(segtype_of_topfile(f), 'ligand')

    def test_model_rule_does_not_catch_real_residue_files(self):
        # `_model` must not match the family files that hold real residues
        self.assertEqual(segtype_of_topfile('toppar_all36_prot_modify_res.str'), 'protein')
        self.assertEqual(segtype_of_topfile('toppar_all36_lipid_prot.str'), 'lipid')
        self.assertEqual(segtype_of_topfile('toppar_all36_carb_glycopeptide.str'), 'glycan')

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


class TestStandaloneMoleculesInPolymerStreams(unittest.TestCase):
    """A residue from a protein or nucleic-acid file that bonds to no neighbour is a molecule, not a
    chain residue.  348 were misclassified this way: dipeptide models, pyridines, CO/O2 for heme,
    coenzymes, free nucleotides."""

    def test_protein_file_standalone_becomes_ligand(self):
        d = derive_segtypes({'CO2': 'toppar_all36_prot_heme_for_new_psf_gen_code_2022.str',
                             'ALA': 'top_all36_prot.rtf'}, standalone={'CO2'})
        self.assertEqual(d, {'ligand': ['CO2'], 'protein': ['ALA']})

    def test_cofactor_streams_give_cofactor(self):
        d = derive_segtypes({'DHF': 'toppar_all36_prot_cofactors.str',
                             'AMP': 'toppar_all36_na_nad_ppi.str',
                             'ADE': 'top_all36_na.rtf'}, standalone={'DHF', 'AMP'})
        self.assertEqual(d, {'cofactor': ['AMP', 'DHF'], 'nucleicacid': ['ADE']})

    def test_only_polymer_families_are_affected(self):
        # sugars and lipids are standalone even when real; the rule must not touch them
        d = derive_segtypes({'BGLC': 'top_all36_carb.rtf', 'POPC': 'top_all36_lipid.rtf'},
                            standalone={'BGLC', 'POPC'})
        self.assertEqual(d, {'glycan': ['BGLC'], 'lipid': ['POPC']})

    def test_without_bond_information_nothing_changes(self):
        d = derive_segtypes({'CO2': 'toppar_all36_prot_heme_for_new_psf_gen_code_2022.str'})
        self.assertEqual(d, {'protein': ['CO2']})
