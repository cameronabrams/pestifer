# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Residue names cut to four characters are restored from what the residue contains.

pestifer used a psfgen alias, BGLC -> BGLCNA, to repair GlcNAc written as BGLC.  psfgen applies an
alias to every residue of that name, so a real beta-glucose (BGLC) was silently built as GlcNAc --
reproduced by rebuilding from example 31's own output.  The stub is ambiguous for about 70 CHARMM
names (BGLCNA, BGLCA, BGLCN all become BGLC), so only the atoms can decide.
"""
import unittest
from pathlib import Path
from types import SimpleNamespace as NS

from pidibble.pdbparse import PDBParser

from pestifer.core.objmanager import ObjManager
from pestifer.molecule.molecule import Molecule  # noqa: F401  -- defines the model Segment refers to
from pestifer.molecule.asymmetricunit import AsymmetricUnit
from pestifer.molecule.chainidmanager import ChainIDManager
from pestifer.molecule.resname_repair import choose_resname, repair_truncated_resnames, _writable_name

GLC = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O2', 'O3', 'O4', 'O5', 'O6'}
GLCNA = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O3', 'O4', 'O5', 'O6', 'N', 'C', 'O', 'CT'}
GLCNA_PDB = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O3', 'O4', 'O5', 'O6', 'N2', 'C7', 'O7', 'C8'}
GLCA = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O2', 'O3', 'O4', 'O5', 'O61', 'O62'}
GLCN = {'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'O3', 'O4', 'O5', 'O6', 'N'}


class TestChooseResname(unittest.TestCase):

    def test_a_real_glucose_keeps_its_name(self):
        self.assertEqual(choose_resname('BGLC', GLC)[0], 'BGLC')

    def test_a_glcnac_is_restored_from_either_atom_naming(self):
        self.assertEqual(choose_resname('BGLC', GLCNA)[0], 'BGLCNA')
        self.assertEqual(choose_resname('BGLC', GLCNA_PDB)[0], 'BGLCNA')

    def test_a_glucuronic_acid_is_restored(self):
        self.assertEqual(choose_resname('BGLC', GLCA)[0], 'BGLCA')

    def test_other_stems(self):
        self.assertEqual(choose_resname('AGAL', GLCNA)[0], 'AGALNA')

    def test_a_name_that_is_not_a_stem_is_untouched(self):
        self.assertEqual(choose_resname('ALA', {'N', 'CA', 'C', 'O', 'CB'}), ('ALA', 'not a truncation stem'))


class TestWritableName(unittest.TestCase):
    """Segment PDBs keep four columns of a residue name, so a long CHARMM name must travel as the PDB
    code aliased to it -- the way a deposit's own sugars reach psfgen."""

    def test_long_names_travel_as_their_pdb_code(self):
        self.assertEqual(_writable_name('BGLCNA'), 'NAG')
        self.assertEqual(_writable_name('BGLCA'), 'BDP')
        self.assertEqual(_writable_name('ANE5AC'), 'SIA')

    def test_short_names_are_written_as_is(self):
        self.assertEqual(_writable_name('BGLC'), 'BGLC')

    def test_a_long_name_with_no_code_cannot_be_written(self):
        self.assertIsNone(_writable_name('BGLCN'))


def _atoms(resname, names, chain='X', resid=1):
    return [NS(resname=resname, name=n, chainID=chain, resid=resid, elem=n[0]) for n in names]


class TestRepairInPlace(unittest.TestCase):

    def test_glucose_and_glcnac_under_the_same_stub_are_told_apart(self):
        atoms = _atoms('BGLC', GLC, resid=1) + _atoms('BGLC', GLCNA, resid=2)
        self.assertEqual(repair_truncated_resnames(atoms), 1)
        self.assertEqual({a.resname for a in atoms if a.resid == 1}, {'BGLC'})
        self.assertEqual({a.resname for a in atoms if a.resid == 2}, {'NAG'})

    def test_an_unwritable_identification_warns_and_leaves_the_name(self):
        atoms = _atoms('BGLC', GLCN)
        with self.assertLogs('pestifer.molecule.resname_repair', level='WARNING') as cm:
            self.assertEqual(repair_truncated_resnames(atoms), 0)
        self.assertEqual({a.resname for a in atoms}, {'BGLC'})
        self.assertTrue(any('BGLCN' in m for m in cm.output), cm.output)


class TestRepairAtIngest(unittest.TestCase):
    """Through AsymmetricUnit, on real coordinates: a GlcNAc cut from example 7's output and a beta-
    glucose from example 31's, both written as BGLC."""

    def test_the_two_come_out_as_different_residues(self):
        pdb = Path(__file__).parent.parent.parent / 'inputs' / '6pti_with_truncated_sugars.pdb'
        pstruct = PDBParser(filepath=pdb, input_format='PDB').parse().parsed
        AU = AsymmetricUnit(parsed=pstruct, chainIDmanager=ChainIDManager(format='PDB'),
                            objmanager=ObjManager(), sourcespecs={'exclude': ["resname == 'HOH'"]},
                            source_format='PDB')
        # glycan residues are renumbered within their segments, so look them up by segtype
        glycans = sorted(r.resname for seg in AU.segments if seg.segtype == 'glycan' for r in seg.residues.data)
        self.assertEqual(glycans, ['BGLC', 'NAG'])
