# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import logging
import os
import unittest

from pestifer.charmmff.charmmffprm import (
    CharmmParamFile,
)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Synthetic test data
# ---------------------------------------------------------------------------

_PRM_TEXT = """\
* Test CHARMM parameter file
*

BONDS
!
CT1  CT2   222.500     1.5380 ! alkane
CT1  NH1   320.000     1.4300 ! peptide
NH1  H     405.000     1.0200 ! NH-H
C    O     620.000     1.2300 ! carbonyl
C    NH1   370.000     1.3450 ! peptide bond

ANGLES
!
CT2  CT1  NH1   67.700   110.000 ! simple angle
H    NH1  C     35.000   120.000 ! amide
O    C    NH1   80.000   122.520 ! carbonyl-amide
CT1  C    O     70.000   121.400 ! backbone
NH1  CT1  CT2   70.000   110.000  22.530   2.179 ! UB angle

DIHEDRALS
!
C    NH1  CT1  CT2      0.2000  1     0.00
C    NH1  CT1  CT2      0.3000  3     0.00
X    C    NH1  X        2.5000  2   180.00 ! wildcard

IMPROPER
!
O    C    NH1  CT1   10.5000  0     0.00

NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
!
CT1    0.000000  -0.020000     2.275000   0.000000  -0.010000     1.900000 ! with 1-4
CT2    0.000000  -0.055000     2.175000 ! no 1-4
NH1    0.000000  -0.200000     1.850000
C      0.000000  -0.110000     2.000000
O      0.000000  -0.120000     1.700000
H      0.000000  -0.046000     0.224500

NBFIX
NH1    O       -0.154919   3.637000 ! test nbfix

END
"""

_STR_TEXT = """\
* Test CHARMM stream file
*

! Topology section (should be ignored by parameter parser)
read rtf card append
* Test topology
*
31 1

MASS  -1  OT   15.999 O
MASS  -1  HT    1.008 H

RESI TIP3  0.000
GROUP
ATOM OH2  OT  -0.834
ATOM H1   HT   0.417
BOND OH2 H1
END

! Parameter section
read para card flex append
*

BONDS
OT  HT   450.000     0.9572 ! water O-H
HT  HT     0.000     1.5139 ! water H-H SHAKE

ANGLES
HT  OT  HT   55.000   104.520 ! water angle

DIHEDRALS
X   OT  HT  X   0.0     1   0.0

NONBONDED nbxmod  5 atom cdiel fshift vatom vdistance vfswitch -
cutnb 14.0 ctofnb 12.0 ctonnb 10.0 eps 1.0 e14fac 1.0 wmin 1.5
!
OT     0.000000  -0.152100     1.768200 ! TIP3P oxygen
HT     0.000000  -0.046000     0.224500 ! TIP3P hydrogen

NBFIX
OT  CT1  -0.120000   3.500000 ! test nbfix from stream

END
"""


class TestCharmmParamFilePrm(unittest.TestCase):
    """Tests for parsing a standalone .prm file (from string)."""

    @classmethod
    def setUpClass(cls):
        cls.p = CharmmParamFile.from_text(_PRM_TEXT)

    def test_bond_count(self):
        self.assertEqual(len(self.p.bonds), 5)

    def test_bond_types(self):
        bt = {(b.type1, b.type2) for b in self.p.bonds}
        self.assertIn(('CT1', 'CT2'), bt)
        self.assertIn(('C', 'O'), bt)

    def test_bond_values(self):
        b = next(b for b in self.p.bonds if b.type1 == 'CT1' and b.type2 == 'CT2')
        self.assertAlmostEqual(b.Kb, 222.5, places=1)
        self.assertAlmostEqual(b.b0, 1.538, places=3)

    def test_angle_count(self):
        self.assertEqual(len(self.p.angles), 5)

    def test_angle_urey_bradley(self):
        ub = next(a for a in self.p.angles if a.Kub is not None)
        self.assertAlmostEqual(ub.Kub, 22.53, places=2)
        self.assertAlmostEqual(ub.S0, 2.179, places=3)

    def test_angle_no_ub(self):
        no_ub = next(a for a in self.p.angles if a.Kub is None)
        self.assertIsNone(no_ub.Kub)
        self.assertIsNone(no_ub.S0)

    def test_dihedral_count(self):
        self.assertEqual(len(self.p.dihedrals), 3)

    def test_dihedral_wildcard(self):
        wc = next(d for d in self.p.dihedrals if d.type1 == 'X')
        self.assertEqual(wc.type4, 'X')
        self.assertAlmostEqual(wc.Kchi, 2.5, places=1)
        self.assertEqual(wc.n, 2)

    def test_improper_count(self):
        self.assertEqual(len(self.p.impropers), 1)

    def test_improper_values(self):
        imp = self.p.impropers[0]
        self.assertEqual(imp.type1, 'O')
        self.assertAlmostEqual(imp.Kpsi, 10.5, places=1)

    def test_nonbonded_count(self):
        self.assertEqual(len(self.p.nonbonded), 6)

    def test_nonbonded_with_14(self):
        nb = self.p.nonbonded['CT1']
        self.assertIsNotNone(nb.ignored14)
        self.assertAlmostEqual(nb.epsilon, -0.020, places=3)
        self.assertAlmostEqual(nb.Rmin_half, 2.275, places=3)
        self.assertAlmostEqual(nb.epsilon14, -0.010, places=3)

    def test_nonbonded_without_14(self):
        nb = self.p.nonbonded['CT2']
        self.assertIsNone(nb.ignored14)

    def test_nbfix_count(self):
        self.assertEqual(len(self.p.nbfix), 1)

    def test_nbfix_values(self):
        n = self.p.nbfix[0]
        self.assertEqual(n.type1, 'NH1')
        self.assertAlmostEqual(n.Emin, -0.154919, places=5)

    def test_no_cmap(self):
        self.assertEqual(len(self.p.cmaps), 0)


class TestCharmmParamFileStr(unittest.TestCase):
    """Tests for parsing parameter sections from a .str file."""

    @classmethod
    def setUpClass(cls):
        cls.s = CharmmParamFile.from_text(_STR_TEXT)

    def test_bond_count(self):
        self.assertEqual(len(self.s.bonds), 2)

    def test_angle_count(self):
        self.assertEqual(len(self.s.angles), 1)

    def test_dihedral_count(self):
        self.assertEqual(len(self.s.dihedrals), 1)

    def test_nonbonded_count(self):
        self.assertEqual(len(self.s.nonbonded), 2)
        self.assertIn('OT', self.s.nonbonded)
        self.assertIn('HT', self.s.nonbonded)

    def test_nbfix_count(self):
        self.assertEqual(len(self.s.nbfix), 1)

    def test_improper_count(self):
        self.assertEqual(len(self.s.impropers), 0)

    def test_no_topology_bleedover(self):
        """Atom types from the topology section (MASS lines, RESI blocks)
        must NOT pollute the parameter records."""
        # OT/HT bonds only come from the BONDS section inside 'read para...END'
        bond_types = {(b.type1, b.type2) for b in self.s.bonds}
        self.assertIn(('OT', 'HT'), bond_types)
        # The topology RESI records must not create spurious parameter entries
        self.assertNotIn('OH2', self.s.nonbonded)


class TestCharmmParamFileMerge(unittest.TestCase):
    """Tests for merging two CharmmParamFile instances."""

    @classmethod
    def setUpClass(cls):
        cls.combined = CharmmParamFile()
        cls.combined.merge(CharmmParamFile.from_text(_PRM_TEXT))
        cls.combined.merge(CharmmParamFile.from_text(_STR_TEXT))

    def test_bond_count(self):
        self.assertEqual(len(self.combined.bonds), 7)  # 5 + 2

    def test_nonbonded_merged(self):
        # PRM has 6 types, STR adds OT and HT
        self.assertIn('OT', self.combined.nonbonded)
        self.assertIn('HT', self.combined.nonbonded)
        self.assertIn('H', self.combined.nonbonded)  # from PRM

    def test_nbfix_count(self):
        self.assertEqual(len(self.combined.nbfix), 2)  # 1 + 1


class TestCharmmParamFileExtract(unittest.TestCase):
    """Tests for extract_for_atomtypes."""

    @classmethod
    def setUpClass(cls):
        combined = CharmmParamFile()
        combined.merge(CharmmParamFile.from_text(_PRM_TEXT))
        combined.merge(CharmmParamFile.from_text(_STR_TEXT))
        cls.atomtypes = {'CT1', 'CT2', 'NH1', 'C', 'O', 'H'}
        cls.minimal = combined.extract_for_atomtypes(cls.atomtypes)

    def test_bonds_filtered(self):
        for b in self.minimal.bonds:
            self.assertIn(b.type1, self.atomtypes)
            self.assertIn(b.type2, self.atomtypes)
        # OT-HT bond must be excluded
        bond_types = {(b.type1, b.type2) for b in self.minimal.bonds}
        self.assertNotIn(('OT', 'HT'), bond_types)

    def test_wildcard_dihedrals_included(self):
        # The X C NH1 X dihedral should survive (X matches anything)
        wc = [d for d in self.minimal.dihedrals if d.type1 == 'X']
        self.assertTrue(len(wc) > 0)

    def test_nonbonded_filtered(self):
        for t in self.minimal.nonbonded:
            self.assertIn(t, self.atomtypes)
        self.assertNotIn('OT', self.minimal.nonbonded)
        self.assertNotIn('HT', self.minimal.nonbonded)

    def test_nbfix_filtered(self):
        # NH1-O nbfix should survive; OT-CT1 should be excluded (OT not in set)
        for n in self.minimal.nbfix:
            self.assertTrue(n.type1 in self.atomtypes or n.type1 == 'X')
            self.assertTrue(n.type2 in self.atomtypes or n.type2 == 'X')
        nbfix_pairs = {(n.type1, n.type2) for n in self.minimal.nbfix}
        self.assertNotIn(('OT', 'CT1'), nbfix_pairs)


class TestCharmmParamFileRoundTrip(unittest.TestCase):
    """Tests for write() / from_file() round-trip."""

    def _round_trip(self, text: str, outfile: str):
        original = CharmmParamFile.from_text(text)
        original.write(outfile)
        reread = CharmmParamFile.from_file(outfile)
        return original, reread

    def test_prm_roundtrip_counts(self):
        original, reread = self._round_trip(_PRM_TEXT, 'test_rt.prm')
        self.assertEqual(len(original.bonds), len(reread.bonds))
        self.assertEqual(len(original.angles), len(reread.angles))
        self.assertEqual(len(original.dihedrals), len(reread.dihedrals))
        self.assertEqual(len(original.impropers), len(reread.impropers))
        self.assertEqual(len(original.nonbonded), len(reread.nonbonded))
        self.assertEqual(len(original.nbfix), len(reread.nbfix))

    def test_prm_roundtrip_bond_values(self):
        original, reread = self._round_trip(_PRM_TEXT, 'test_rt2.prm')
        for b1, b2 in zip(
            sorted(original.bonds, key=lambda b: (b.type1, b.type2)),
            sorted(reread.bonds, key=lambda b: (b.type1, b.type2)),
        ):
            self.assertEqual(b1.type1, b2.type1)
            self.assertEqual(b1.type2, b2.type2)
            self.assertAlmostEqual(b1.Kb, b2.Kb, places=2)
            self.assertAlmostEqual(b1.b0, b2.b0, places=3)

    def test_prm_roundtrip_nonbonded_values(self):
        original, reread = self._round_trip(_PRM_TEXT, 'test_rt3.prm')
        for t in original.nonbonded:
            self.assertIn(t, reread.nonbonded)
            self.assertAlmostEqual(
                original.nonbonded[t].epsilon, reread.nonbonded[t].epsilon, places=5
            )
            self.assertAlmostEqual(
                original.nonbonded[t].Rmin_half, reread.nonbonded[t].Rmin_half, places=5
            )

    def test_str_roundtrip_counts(self):
        """Ensure round-trip from .str source produces a valid .prm."""
        original, reread = self._round_trip(_STR_TEXT, 'test_str_rt.prm')
        self.assertEqual(len(original.bonds), len(reread.bonds))
        self.assertEqual(len(original.nonbonded), len(reread.nonbonded))
