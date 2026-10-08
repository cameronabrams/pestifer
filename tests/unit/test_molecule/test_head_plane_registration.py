# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A leaflet is registered at the polar/apolar interface, not at each species' head marker.

Where a species' head group sits in a leaflet is fixed by two things, and only one of them is a
property of the conformer ensemble:

``iface_off``
    anchor -> the molecule's lowest polar (N/O/P) atom.  Rigid headgroup geometry: a phosphate
    sits ~5.2 A above POPC's deepest ester oxygen, and 0.00 A above CHL1's, whose anchor *is* its
    3-hydroxyl.  This sets the RELATIVE heights of two species.
``depth``
    interface -> mean tail z.  How far the tails actually reach, which legitimately depends on
    how melted the ensemble is.  Composition-weighted, it sets the ABSOLUTE thickness.

The old rule was one number, ``anchor -> mean tail z``, taken per species -- so the relative
placement of two species depended on their ensembles' extension, and a melted phospholipid is not
comparable to a rigid sterol that way.  Measured 2026-10-08: with an honest POPC ensemble, POPC's
offset came out 15.76 A against CHL1's 16.93, putting the sterol hydroxyl 1.17 A ABOVE the
phosphate plane.

It had only ever looked right because of a defect elsewhere.  The shipped POPC ensemble was seven
copies of one extended draw -- the ``run_mc`` sample checkpoint sat inside the accepted-move
branch -- which biased its mean 6.6 A high and accidentally cleared the sterol.  Fixing the
sampler exposed the registration.  Two wrongs; the collection was the one that got fixed first.

The measured consequence, 50:50 POPC:CHL1: built P-P 46.75 -> 39.88 A, against ~38 A for POPC,
and the sterol hydroxyl at -5.19 A, inside the -4.4 to -5.1 A measured on equilibrated ex17
patches.
"""
import unittest

import numpy as np

from pestifer.molecule.bilayer import _polar_rows

_PDB = ('ATOM  {:5d} {:<4s} RES A   1    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00      SEG')


def _line(i, name):
    return _PDB.format(i, name, 0.0, 0.0, 0.0)


class TestPolarRows(unittest.TestCase):
    """The interface is found by element, read off the atom NAME -- the conformer PDBs leave the
    element column blank.  Measured across all 277 shipped entries and 2698 conformers: the first
    character of every atom name is one of C, H, O, P, N -- no leading digits, no two-letter
    metals -- so a first-character test is exact for this collection."""

    def test_it_finds_N_O_and_P_and_nothing_else(self):
        names = ['N', 'C12', 'H12A', 'P', 'O13', 'C21', 'O22', 'C318', 'HX']
        rows = _polar_rows([_line(i + 1, n) for i, n in enumerate(names)])
        self.assertEqual([names[i] for i in rows], ['N', 'P', 'O13', 'O22'])

    def test_a_cardiolipin_CA16_is_carbon_not_calcium(self):
        """The trap this file is one step away from: a reader guessing an element harder than
        this would read CA16 as calcium. Pinned because the names are all there is to read."""
        rows = _polar_rows([_line(1, 'CA16'), _line(2, 'CA9'), _line(3, 'OA1')])
        self.assertEqual(rows, [2])

    def test_no_polar_atoms_is_empty_not_an_error(self):
        self.assertEqual(_polar_rows([_line(1, 'C1'), _line(2, 'H1')]), [])


class TestTheRegistrationIsEnsembleIndependent(unittest.TestCase):
    """The property the fix exists to establish, stated directly.

    ``iface_off`` must be a function of the molecule, not of the draw.  If it ever starts
    depending on the ensemble, two species' relative placement starts moving whenever the
    collection is regenerated -- which is exactly how a sterol ended up above the phosphates.
    """

    @staticmethod
    def _conformer(anchor_z, iface_z, tail_z, jitter=0.0):
        """A toy lipid: P at anchor_z, an ester O at iface_z, a tail carbon at tail_z."""
        names = ['P', 'O21', 'C218']
        lines = [_line(i + 1, n) for i, n in enumerate(names)]
        coords = np.array([[0.0, 0.0, anchor_z],
                           [0.0, 0.0, iface_z],
                           [0.0, 0.0, tail_z + jitter]])
        return coords, lines

    def test_iface_off_does_not_move_when_the_tails_melt(self):
        """Same headgroup, wildly different tail extension -> same iface_off."""
        offs = []
        for tail_z in (-20.0, -14.0, -8.0):      # extended -> melted
            coords, lines = self._conformer(anchor_z=5.0, iface_z=0.0, tail_z=tail_z)
            offs.append(float(coords[0, 2]) - float(coords[_polar_rows(lines), 2].min()))
        self.assertEqual(len(set(offs)), 1, f'iface_off moved with tail extension: {offs}')
        self.assertAlmostEqual(offs[0], 5.0)

    def test_the_OLD_statistic_DOES_move_when_the_tails_melt(self):
        """NEGATIVE CONTROL: the quantity that was being used is ensemble-dependent.

        Without this, the test above could pass against any rigid-ish toy molecule and say
        nothing about why the statistic was changed.
        """
        offs = []
        for tail_z in (-20.0, -14.0, -8.0):
            coords, _ = self._conformer(anchor_z=5.0, iface_z=0.0, tail_z=tail_z)
            offs.append(float(coords[0, 2]) - float(coords[2, 2]))   # anchor -> tail, the old rule
        self.assertEqual(len(set(offs)), 3, 'the old statistic was supposed to track the tails')

    def test_a_rigid_sterol_lands_BELOW_a_phospholipid_head(self):
        """The sign that was wrong. CHL1's anchor is its own hydroxyl, so its iface_off is 0;
        any phospholipid's phosphate sits several A above its deepest ester oxygen."""
        popc, popc_lines = self._conformer(anchor_z=5.19, iface_z=0.0, tail_z=-15.5)
        chl1, chl1_lines = self._conformer(anchor_z=0.0, iface_z=0.0, tail_z=-16.93)
        off_popc = float(popc[0, 2]) - float(popc[_polar_rows(popc_lines), 2].min())
        off_chl1 = float(chl1[0, 2]) - float(chl1[_polar_rows(chl1_lines), 2].min())
        self.assertLess(off_chl1 - off_popc, -1.0,
                        'the sterol hydroxyl is not below the phosphate plane')
        self.assertAlmostEqual(off_chl1 - off_popc, -5.19, places=2)


if __name__ == '__main__':
    unittest.main()
