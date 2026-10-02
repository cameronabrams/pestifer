# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Tests for :func:`pestifer.util.density_convergence.membrane_pp_thickness`.

Phosphate-to-phosphate thickness is the structural observable the membrane-equilibrate gate was
missing: area overshoots on release (50.00 -> 54.53 -> 47.71 on ex17 patchB) and its drift passes
through zero at the turn while the bilayer is still ordering, so an area-only test certifies a
structure that has not finished moving.  Thickness does not turn (2 direction reversals against
area's 5) and tracks chain order almost exactly (corr 0.955), which is why it is worth measuring.

Everything here is synthetic and exact: a bilayer is written with its phosphorus atoms at known z,
so the answer is known to the digit rather than approximated.  No NAMD, no real structure.
"""
import os
import tempfile
import unittest

from pestifer.util.density_convergence import (MembraneThickness, membrane_leaflet_geometry,
                                               membrane_pp_thickness)

P_MASS = 30.974
C_MASS = 12.011
N_MASS = 14.007


def _write_system(d, molecules):
    """Write a (psf, pdb) pair.  ``molecules`` is a list of ``(resname, [(name, mass, z), ...])``."""
    psf, pdb = os.path.join(d, 's.psf'), os.path.join(d, 's.pdb')
    atoms = []
    for mol_i, (resname, spec) in enumerate(molecules, start=1):
        for name, mass, z in spec:
            atoms.append((f'L{mol_i % 2}', mol_i, resname, name, mass, z))
    with open(psf, 'w') as f:
        f.write(f'PSF EXT\n\n{len(atoms):8d} !NATOM\n')
        for i, (seg, resid, resn, name, mass, _z) in enumerate(atoms, start=1):
            f.write(f'{i:8d} {seg:8s} {resid:<8d} {resn:8s} {name:8s} {name:8s} '
                    f'0.00 {mass:10.4f}  0\n')
    with open(pdb, 'w') as f:
        for i, (seg, resid, resn, name, _m, z) in enumerate(atoms, start=1):
            f.write(f'ATOM  {i:5d} {name:<4s} {resn:<4s}{resid:5d}    '
                    f'{0.0:8.3f}{0.0:8.3f}{z:8.3f}  1.00  0.00\n')
        f.write('END\n')
    return psf, pdb


def _phospholipid(resname, p_z, sign):
    """A minimal phospholipid: one phosphorus at ``p_z`` plus two tail carbons toward the midplane."""
    return (resname, [('P', P_MASS, p_z),
                      ('C21', C_MASS, p_z - sign * 4.0),
                      ('C22', C_MASS, p_z - sign * 8.0)])


def _bilayer(p_upper=20.0, p_lower=-20.0, n_each=4, resname='POPC'):
    return ([_phospholipid(resname, p_upper, +1) for _ in range(n_each)]
            + [_phospholipid(resname, p_lower, -1) for _ in range(n_each)])


class TestPPThickness(unittest.TestCase):

    def _measure(self, molecules):
        with tempfile.TemporaryDirectory() as d:
            psf, pdb = _write_system(d, molecules)
            return membrane_pp_thickness(psf, pdb)

    def test_a_clean_bilayer_measures_exactly(self):
        r = self._measure(_bilayer(20.0, -20.0))
        self.assertAlmostEqual(r.thickness, 40.0, places=6)
        self.assertAlmostEqual(r.z_upper, 20.0, places=6)
        self.assertAlmostEqual(r.z_lower, -20.0, places=6)
        self.assertEqual((r.n_upper, r.n_lower), (4, 4))
        self.assertEqual((r.n_no_phosphorus, r.n_unresolved), (0, 0))

    def test_an_asymmetric_bilayer_is_not_assumed_symmetric(self):
        """The midplane is mass-weighted, not zero; thickness must not be 2*|z_upper|."""
        r = self._measure(_bilayer(30.0, -10.0))
        self.assertAlmostEqual(r.thickness, 40.0, places=6)

    def test_thickness_is_a_difference_of_leaflet_means(self):
        mols = ([_phospholipid('POPC', 19.0, +1), _phospholipid('POPC', 21.0, +1)]
                + [_phospholipid('POPC', -20.0, -1), _phospholipid('POPC', -20.0, -1)])
        r = self._measure(mols)
        self.assertAlmostEqual(r.thickness, 40.0, places=6)

    def test_cholesterol_is_excluded_not_an_error(self):
        """CHL1 is 43-47% of ex17's leaflets and has no phosphorus.  P-P is a phospholipid
        measure, so sterols are counted and skipped -- never guessed at."""
        chol = ('CHL1', [('C1', C_MASS, 15.0), ('O1', 15.999, 18.0)])
        r = self._measure(_bilayer() + [chol, chol])
        self.assertAlmostEqual(r.thickness, 40.0, places=6)
        self.assertEqual(r.n_no_phosphorus, 2)
        self.assertEqual((r.n_upper, r.n_lower), (4, 4), 'sterols must not enter the P counts')

    def test_a_leaflet_with_no_phospholipid_is_undefined_not_zero(self):
        """NEGATIVE CONTROL: the observable must refuse to answer rather than invent a number.

        A caller that gates on thickness has to treat this as "do not certify".  Returning 0.0,
        or silently measuring from one leaflet, is the failure this test exists to prevent.
        """
        sterol_only = ('CHL1', [('C1', C_MASS, 20.0)])
        mols = [_phospholipid('POPC', -20.0, -1) for _ in range(3)] + [sterol_only]
        r = self._measure(mols)
        self.assertIsNone(r.thickness)
        self.assertIn('upper', r.note)
        self.assertEqual(r.n_lower, 3)

    def test_both_leaflets_empty_is_reported_as_such(self):
        r = self._measure([('CHL1', [('C1', C_MASS, 10.0)]), ('CHL1', [('C1', C_MASS, -10.0)])])
        self.assertIsNone(r.thickness)
        self.assertIn('both leaflets', r.note)

    def test_a_multi_phosphorus_lipid_resolves_to_the_headgroup_P(self):
        """PIP2 carries three phosphorus atoms; averaging them would bias the headgroup plane."""
        pip2_up = ('PIP2', [('P', P_MASS, 20.0), ('P4', P_MASS, 26.0), ('P5', P_MASS, 28.0),
                            ('C21', C_MASS, 16.0)])
        r = self._measure(_bilayer() + [pip2_up])
        self.assertAlmostEqual(r.thickness, 40.0, places=6,
                               msg='P4/P5 leaked into the headgroup mean')
        self.assertEqual(r.n_upper, 5)
        self.assertEqual(r.n_unresolved, 0)

    def test_an_ambiguous_multi_phosphorus_lipid_is_excluded_and_counted(self):
        """Cardiolipin's two phosphates are equivalent -- neither is named `P`.  Excluding it is
        honest; averaging them into the plane would not be, and would never be visible."""
        cardio = ('TOCL', [('P1', P_MASS, 24.0), ('P2', P_MASS, 28.0), ('C21', C_MASS, 16.0)])
        r = self._measure(_bilayer() + [cardio])
        self.assertAlmostEqual(r.thickness, 40.0, places=6)
        self.assertEqual(r.n_unresolved, 1)
        self.assertEqual(r.n_upper, 4, 'the ambiguous lipid must not be counted as resolved')

    def test_phosphorus_is_found_by_mass_not_by_name(self):
        """No per-lipid atom-name table: a phosphorus named anything still counts, and a
        carbon named `P` does not."""
        odd = ('XYZ', [('PA', P_MASS, 20.0), ('C1', C_MASS, 16.0)])
        r = self._measure([odd] + [_phospholipid('POPC', -20.0, -1)])
        self.assertAlmostEqual(r.thickness, 40.0, places=6)

        decoy = ('XYZ', [('P', C_MASS, 20.0), ('C1', C_MASS, 16.0)])
        r2 = self._measure([decoy] + [_phospholipid('POPC', -20.0, -1)])
        self.assertIsNone(r2.thickness, 'a carbon named P must not be read as phosphorus')

    def test_a_dipping_headgroup_stays_with_its_own_leaflet(self):
        """The claim `_lipid_molecules` makes, which nothing else pinned.

        A lipid's phosphate can dip past the midplane while the body of the molecule is plainly
        in the upper leaflet.  Assigning by the *molecule's* mass-weighted mean keeps it upper --
        correct, since it is an upper-leaflet lipid with a transiently low headgroup, and its P
        belongs in the upper mean.  Assigning by the phosphorus' own z would move it to the lower
        leaflet and corrupt both means at once.

        Caught by mutation testing: swapping the assignment to use the P's z left all twelve
        other tests green, because every synthetic lipid here has its P at the extreme.
        """
        dipping = ('POPC', [('P', P_MASS, -1.0), ('C21', C_MASS, 12.0), ('C22', C_MASS, 16.0)])
        r = self._measure(_bilayer(20.0, -20.0, n_each=4) + [dipping])
        self.assertEqual(r.n_upper, 5, 'the dipping lipid was assigned by its P, not its body')
        self.assertEqual(r.n_lower, 4)
        self.assertAlmostEqual(r.z_upper, (20.0 * 4 + -1.0) / 5, places=6)

    def test_not_a_membrane_raises(self):
        from pestifer.util.densityprofile import WATER_RESNAMES
        water = (sorted(WATER_RESNAMES)[0], [('OH2', 15.999, 0.0)])
        with self.assertRaises(ValueError):
            self._measure([water, water])


class TestLeafletGeometryStillWorks(unittest.TestCase):
    """Regression on the refactor: `membrane_leaflet_geometry` now shares `_lipid_molecules` with
    the thickness measurement, so its own counting must be unchanged."""

    def test_lipid_counts_survive_the_shared_helper(self):
        with tempfile.TemporaryDirectory() as d:
            psf, pdb = _write_system(d, _bilayer(20.0, -20.0, n_each=5))
            g = membrane_leaflet_geometry(psf, pdb)
        self.assertEqual((g.n_upper, g.n_lower), (5, 5))
        self.assertEqual(g.n_lipids_total, 10)
        self.assertAlmostEqual(g.midplane_z, 0.0, places=6)

    def test_apl_is_still_protein_corrected_per_leaflet(self):
        with tempfile.TemporaryDirectory() as d:
            psf, pdb = _write_system(d, _bilayer(20.0, -20.0, n_each=4))
            g = membrane_leaflet_geometry(psf, pdb)
        self.assertAlmostEqual(g.apl_lower(400.0), 100.0, places=6)
        self.assertAlmostEqual(g.apl_mean(400.0), 100.0, places=6)


if __name__ == '__main__':
    unittest.main()
