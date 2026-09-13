# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Small-molecule model compounds are ligands, whichever force-field stream defines them.

pestifer derives a residue's segtype from the FILE that defines it.  CHARMM keeps model compounds
-- small standalone molecules used to parameterize a family -- in toppar_all36_{prot,carb,lipid,
na}_model.str, so 193 of them derived as protein, glycan, lipid or nucleic acid.  Among them were
acetate (example 5's crystallographic ligand), acetone (example 24's solvent, 1,753 molecules
classified as sugar), and the methanol and ethanol pestifer ships as solvent boxes, both classified
as protein.  A `_model` rule now derives them as ligand, matching DMSO and acetonitrile.  The
exceptions are curated: MLYS is a real chain residue, and three phosphates stay lipid.
"""
import unittest

from pestifer.core.labels import Labels, _load_derived_segtypes


class TestModelCompoundsAreLigands(unittest.TestCase):

    def test_named_cases_are_ligands(self):
        cases = {'ACET': 'acetate, example 5', 'ACO': 'acetone, example 24',
                 'MEOH': 'methanol, shipped solvent box', 'ETOH': 'ethanol, shipped solvent box',
                 'PRO2': '2-propanol', 'EGLY': 'ethylene glycol', 'NUTA': 'model nucleotide'}
        for r, desc in cases.items():
            with self.subTest(residue=r):
                self.assertEqual(Labels.segtype_of_resname.get(r), 'ligand', f'{r} ({desc})')

    def test_they_come_from_the_derivation_not_hand_curation(self):
        # the rule, not a list of names, is what classifies them -- so the next model compound
        # CHARMM adds is covered too.  The curated set is read in a FRESH interpreter: the curated
        # lists are mutable at runtime (a config's psfgen.segtypes extends them, and example 5's
        # asks for `other: [ACET, ACT]`), so in-process it depends on which tests ran first.
        import json, subprocess, sys
        derived_ligands = set(_load_derived_segtypes().get('ligand', []))
        curated = set(json.loads(subprocess.run(
            [sys.executable, '-c', 'import json; from pestifer.core.labels import curated_resname_set; '
                                   'print(json.dumps(sorted(curated_resname_set())))'],
            capture_output=True, text=True, check=True).stdout))
        self.assertIn('HEM', curated)   # the subprocess really read the curated table
        for r in ('ACET', 'ACO', 'MEOH', 'ETOH'):
            self.assertIn(r, derived_ligands)
            self.assertNotIn(r, curated)

    def test_mlys_is_a_chain_residue_and_stays_protein(self):
        # defined in prot_model, but it bonds -C/+N to its neighbors
        self.assertEqual(Labels.segtype_of_resname.get('MLYS'), 'protein')
        self.assertIn('MLYS', Labels.segtypes['protein']['resnames'])

    def test_curated_phosphates_stay_lipid(self):
        for r in ('DMP', 'MP_1', 'MP_2'):
            self.assertEqual(Labels.segtype_of_resname.get(r), 'lipid')

    def test_real_lipids_and_amino_acids_are_untouched(self):
        for r, st in (('CHL1', 'lipid'), ('POPC', 'lipid'), ('CYSP', 'protein'), ('SEP', 'protein'),
                      ('DMSO', 'ligand'), ('ACN', 'ligand')):
            self.assertEqual(Labels.segtype_of_resname.get(r), st, r)


if __name__ == '__main__':
    unittest.main()
