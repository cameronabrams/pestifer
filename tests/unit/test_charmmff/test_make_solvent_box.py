import os
import tempfile
import unittest

import numpy as np

from pestifer.charmmff.make_solvent_box import (
    box_edge_for_density, cubic_lattice_sites, pack_cubic,
    write_box_pdb, atom_lines_of, coords_of,
)


class TestSolventBoxGeometry(unittest.TestCase):
    def test_edge_matches_known_water_box(self):
        # 216 TIP3P waters at 1.0 g/cc -> ~18.6 A (VMD's built-in 216-water box is ~18.77)
        e = box_edge_for_density(216, 18.015, 1.0)
        self.assertAlmostEqual(e, 18.63, places=1)

    def test_edge_scaling_and_validation(self):
        # 8x the molecules at the same density -> 2x the edge
        e1 = box_edge_for_density(216, 18.0, 1.0)
        e8 = box_edge_for_density(216 * 8, 18.0, 1.0)
        self.assertAlmostEqual(e8 / e1, 2.0, places=6)
        for bad in [(0, 18, 1), (216, 0, 1), (216, 18, 0)]:
            with self.assertRaises(ValueError):
                box_edge_for_density(*bad)

    def test_lattice_sites_fill_and_bound(self):
        e = 20.0
        s = cubic_lattice_sites(216, e)         # perfect cube -> 6x6x6
        self.assertEqual(s.shape, (216, 3))
        self.assertTrue((s >= 0).all() and (s < e).all())
        self.assertEqual(cubic_lattice_sites(100, e).shape, (100, 3))   # non-cube -> next cube, truncated
        self.assertEqual(cubic_lattice_sites(1, e).shape, (1, 3))

    def test_pack_is_rigid_and_deterministic(self):
        mol = np.array([[0.0, 0, 0], [1.5, 0, 0], [0, 1.0, 0]])
        e = box_edge_for_density(64, 40.0, 0.8)
        p = pack_cubic(mol, 64, e, seed=7)
        self.assertEqual(p.shape, (64, 3, 3))
        # random rotation preserves internal geometry for every copy
        for i in range(64):
            self.assertAlmostEqual(np.linalg.norm(p[i, 1] - p[i, 0]), 1.5, places=5)
            self.assertAlmostEqual(np.linalg.norm(p[i, 2] - p[i, 0]), 1.0, places=5)
        self.assertTrue(np.allclose(p, pack_cubic(mol, 64, e, seed=7)))   # deterministic

    def test_write_box_pdb_roundtrips(self):
        template = [
            "ATOM      1  O   GRL X   1       0.000   0.000   0.000  1.00  0.00      X    O\n",
            "ATOM      2  C1  GRL X   1       1.400   0.000   0.000  1.00  0.00      X    C\n",
        ]
        mol = coords_of(template)
        self.assertEqual(mol.shape, (2, 3))
        placed = pack_cubic(mol, 8, box_edge_for_density(8, 32.0, 0.79), seed=0)
        d = tempfile.mkdtemp()
        out = os.path.join(d, 'box.pdb')
        segids = write_box_pdb(placed, template, out)
        self.assertEqual(len(segids), 8)
        self.assertEqual(len(set(segids)), 8)                 # unique per molecule
        lines = atom_lines_of(out)
        self.assertEqual(len(lines), 16)                       # 8 molecules x 2 atoms
        self.assertEqual(lines[0][72:76].strip(), segids[0])   # segid written in the right column
        self.assertEqual(lines[2][72:76].strip(), segids[1])
        # coordinates round-trip back to the placed positions
        self.assertTrue(np.allclose(coords_of(lines[:2]), placed[0], atol=1e-3))


if __name__ == '__main__':
    unittest.main()
