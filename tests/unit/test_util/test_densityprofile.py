"""Unit tests for pestifer.util.densityprofile."""
import os
import tempfile
import unittest

import numpy as np

from pestifer.util.densityprofile import (
    DensityProfile, classify_species, AMU_PER_A3_TO_G_PER_CC)


def _pdb_line(serial, name, resname, resid, x, y, z):
    """Build a fixed-column PDB ATOM record (z occupies columns 47-54)."""
    line = [' '] * 80

    def put(s, start):  # start is 1-indexed
        for i, ch in enumerate(s):
            line[start - 1 + i] = ch

    put('ATOM', 1)
    put(f'{serial:>5}', 7)
    put(f'{name:<4}', 13)
    put(f'{resname:<4}', 18)
    put(f'{resid:>4}', 23)
    put(f'{x:8.3f}', 31)
    put(f'{y:8.3f}', 39)
    put(f'{z:8.3f}', 47)
    return ''.join(line).rstrip() + '\n'


def _write_system(dirpath, atoms, a_x=10.0, b_y=10.0, c_z=100.0):
    """atoms: list of (resname, name, mass, z).  Writes base.psf/.pdb/.xsc."""
    base = os.path.join(dirpath, 'sys')
    n = len(atoms)
    with open(base + '.psf', 'w') as f:
        f.write('PSF\n\n')
        f.write(f'{n:>8} !NATOM\n')
        for i, (resname, name, mass, z) in enumerate(atoms, start=1):
            seg = 'X'
            f.write(f'{i:>8} {seg:<4} {1:<4} {resname:<6} {name:<6} '
                    f'{name:<6} {0.0:>10.6f} {mass:>10.4f} {0:>5}\n')
    with open(base + '.pdb', 'w') as f:
        for i, (resname, name, mass, z) in enumerate(atoms, start=1):
            f.write(_pdb_line(i, name, resname, 1, 0.0, 0.0, z))
        f.write('END\n')
    with open(base + '.xsc', 'w') as f:
        f.write('# NAMD extended system configuration\n')
        f.write(f'0 {a_x} 0 0 0 {b_y} 0 0 0 {c_z} 0 0 0\n')
    return base


class TestClassify(unittest.TestCase):
    def test_classify(self):
        res = np.array(['TIP3', 'POT', 'CLA', 'ALA', 'DMPC', 'CHL1'], dtype=object)
        cls = classify_species(res)
        self.assertEqual(list(cls),
                         ['water', 'ion', 'ion', 'protein', 'lipid', 'lipid'])


class TestDensityProfile(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        # midplane of the lipids will be at z=50; cell c_z=100 -> edges [0,100] about it
        self.atoms = [
            ('DMPC', 'C1', 12.0, 48.0),
            ('DMPC', 'C2', 12.0, 52.0),
            ('CHL1', 'C1', 12.0, 46.0),
            ('CHL1', 'C2', 12.0, 54.0),
            ('TIP3', 'OH2', 16.0, 20.0),
            ('TIP3', 'OH2', 16.0, 80.0),
            ('ALA', 'CA', 12.0, 50.0),
            ('POT', 'POT', 39.0, 10.0),
        ]
        self.base = _write_system(self.tmp, self.atoms, a_x=10.0, b_y=10.0, c_z=100.0)

    def test_lipid_resnames(self):
        dp = DensityProfile(self.base + '.psf', self.base + '.pdb', self.base + '.xsc')
        self.assertEqual(dp.lipid_resnames, ['CHL1', 'DMPC'])
        self.assertAlmostEqual(dp.area, 100.0)
        self.assertAlmostEqual(dp.c_z, 100.0)

    def test_mass_conservation(self):
        """Integrating each species' density recovers its total mass."""
        dp = DensityProfile(self.base + '.psf', self.base + '.pdb', self.base + '.xsc')
        centers, profiles = dp.compute(dz=1.0)
        slab_vol = dp.area * (dp.c_z / len(centers))
        # total water mass = 2 * 16
        recovered = profiles['water'].sum() * slab_vol / AMU_PER_A3_TO_G_PER_CC
        self.assertAlmostEqual(recovered, 32.0, places=5)
        # total lipid mass = 4 * 12
        recovered = profiles['lipid'].sum() * slab_vol / AMU_PER_A3_TO_G_PER_CC
        self.assertAlmostEqual(recovered, 48.0, places=5)

    def test_pbc_wrapping(self):
        """An atom beyond the box edge wraps in (mass is not dropped)."""
        atoms = list(self.atoms)
        # place an extra water at z = midplane + c_z/2 + 5 = 50+50+5 = 105 (past the top)
        atoms.append(('TIP3', 'OH2', 16.0, 105.0))
        base = _write_system(self.tmp, atoms, a_x=10.0, b_y=10.0, c_z=100.0)
        dp = DensityProfile(base + '.psf', base + '.pdb', base + '.xsc')
        centers, profiles = dp.compute(dz=1.0)
        slab_vol = dp.area * (dp.c_z / len(centers))
        recovered = profiles['water'].sum() * slab_vol / AMU_PER_A3_TO_G_PER_CC
        # 3 waters * 16 -- the wrapped one is retained
        self.assertAlmostEqual(recovered, 48.0, places=5)

    def test_components_sum_to_total(self):
        dp = DensityProfile(self.base + '.psf', self.base + '.pdb', self.base + '.xsc')
        centers, profiles = dp.compute(dz=1.0, lipid_components=True)
        self.assertIn('lipid:DMPC', profiles)
        self.assertIn('lipid:CHL1', profiles)
        comp_sum = profiles['lipid:DMPC'] + profiles['lipid:CHL1']
        np.testing.assert_allclose(comp_sum, profiles['lipid'], atol=1e-9)

    def test_namdbin_coor(self):
        """A NAMD binary .coor frame yields the same profile as the PDB."""
        import struct
        z = np.array([a[3] for a in self.atoms], dtype='<f8')
        xyz = np.zeros((len(z), 3), dtype='<f8')
        xyz[:, 2] = z
        with open(self.base + '.coor', 'wb') as f:
            f.write(struct.pack('<i', len(z)))
            f.write(xyz.tobytes())
        dp_pdb = DensityProfile(self.base + '.psf', self.base + '.pdb', self.base + '.xsc')
        dp_bin = DensityProfile(self.base + '.psf', self.base + '.coor', self.base + '.xsc')
        _, p_pdb = dp_pdb.compute(dz=1.0)
        _, p_bin = dp_bin.compute(dz=1.0)
        np.testing.assert_allclose(p_bin['lipid'], p_pdb['lipid'], atol=1e-9)


if __name__ == '__main__':
    unittest.main()
