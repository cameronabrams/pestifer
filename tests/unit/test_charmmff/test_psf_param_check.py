# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the incoming-PSF force-field-consistency checker."""
import unittest
from pathlib import Path
import tempfile

from pestifer.charmmff.charmmffprm import CharmmParamFile
from pestifer.charmmff.psf_param_check import check_psf_parameters, format_missing
from pestifer.psfutil.psfcontents import PSFContents


# A synthetic parameter set covering exactly the terms the crafted PSF below uses.
# The dihedral is a wildcard 'X CT1 C X' (tests central-pair wildcard matching); the
# improper is a wildcard 'C X X O' (tests permissive improper wildcard matching).
_PRM_TEXT = """\
* synthetic test parameters
*

BONDS
NH1  CT1   300.0   1.45
CT1  C     250.0   1.49
C    O     620.0   1.23

ANGLES
NH1  CT1  C    50.0   110.0
CT1  C    O    80.0   121.0

DIHEDRALS
X    CT1  C    X     0.2000  1   0.00

IMPROPER
C    X    X    O    120.0   0   0.00

NONBONDED nbxmod 5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 16.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
NH1   0.0  -0.20   1.85
CT1   0.0  -0.032  2.00
C     0.0  -0.11   2.00
O     0.0  -0.12   1.70
"""

# A minimal 4-atom alanine-like residue: one dihedral (NH1-CT1-C-O) and one improper.
_PSF_TEMPLATE = """\
PSF EXT

       3 !NTITLE
 REMARKS test system
 REMARKS topology top_all36_prot.rtf
 REMARKS segment PROA { first NTER; last CTER; auto angles dihedrals }

       4 !NATOM
       1 PROA     1        RESNM     N        NH1   -0.470000       14.0070           0
       2 PROA     1        RESNM     CA       CT1    0.070000       12.0110           0
       3 PROA     1        RESNM     C        C      0.510000       12.0110           0
       4 PROA     1        RESNM     O        O     -0.510000       15.9990           0

       3 !NBOND: bonds
       1       2       2       3       3       4

       2 !NTHETA: angles
       1       2       3       2       3       4

       1 !NPHI: dihedrals
       1       2       3       4

       1 !NIMPHI: impropers
       3       1       2       4

       0 !NDON: donors

       0 !NACC: acceptors

       0 !NNB
"""


def _write_psf(dirpath, resname='ALA'):
    p = Path(dirpath) / f'{resname}_test.psf'
    p.write_text(_PSF_TEMPLATE.replace('RESNM', resname))
    return str(p)


class TestPsfParamCheck(unittest.TestCase):

    def setUp(self):
        self.param = CharmmParamFile.from_text(_PRM_TEXT)
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_all_terms_resolve(self):
        psf = PSFContents(_write_psf(self.dir))
        missing = check_psf_parameters(psf, self.param)
        self.assertFalse(missing.any())
        # in particular the wildcard dihedral and improper must be considered resolved
        self.assertEqual(missing.dihedrals, [])
        self.assertEqual(missing.impropers, [])

    def test_missing_atom_type_is_named_with_residue(self):
        psf = PSFContents(_write_psf(self.dir))
        self.param.nonbonded.pop('CT1')
        missing = check_psf_parameters(psf, self.param)
        self.assertTrue(missing.any())
        self.assertEqual(len(missing.atomtypes), 1)
        term, where = missing.atomtypes[0]
        self.assertEqual(term, 'CT1')
        self.assertEqual(where, 'ALA PROA1')
        msg = format_missing(missing, 'feb26')
        self.assertIn("'CT1'", msg)
        self.assertIn('ALA PROA1', msg)
        self.assertIn('feb26', msg)

    def test_missing_bond_is_flagged(self):
        psf = PSFContents(_write_psf(self.dir))
        self.param.bonds = [b for b in self.param.bonds
                            if tuple(sorted((b.type1, b.type2))) != ('CT1', 'NH1')]
        missing = check_psf_parameters(psf, self.param)
        self.assertTrue(any('NH1' in t and 'CT1' in t for t, _ in missing.bonds))

    def test_missing_wildcard_dihedral_is_flagged(self):
        psf = PSFContents(_write_psf(self.dir))
        self.param.dihedrals = []   # remove the only (wildcard) dihedral
        missing = check_psf_parameters(psf, self.param)
        self.assertEqual(len(missing.dihedrals), 1)

    def test_unknown_resname_parses_without_crashing(self):
        # A resname unknown to pestifer's segtype table must not crash PSF parsing;
        # its topology is still readable and its terms still checkable.
        psf = PSFContents(_write_psf(self.dir, resname='LIG'))
        self.assertEqual(len(psf.atoms), 4)
        self.assertEqual(psf.atoms.data[0].segtype, '')   # classifies blank, not a crash
        missing = check_psf_parameters(psf, self.param)
        self.assertFalse(missing.any())


if __name__ == '__main__':
    unittest.main()
