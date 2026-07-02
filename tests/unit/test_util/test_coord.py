import os
import tempfile
import numpy as np
import unittest

from pestifer.util.coord import rotate_points_about_axis, pdb_replace_coords


class TestRotatePointsAboutAxis(unittest.TestCase):
    def test_180_about_z_through_origin(self):
        pts = np.array([[1.0, 0.0, 5.0], [0.0, 2.0, -3.0]])
        r = rotate_points_about_axis(pts, pivot=np.zeros(3), axis=[0, 0, 1], degrees=180)
        # a 180 deg rotation about z negates x and y, leaves z
        self.assertTrue(np.allclose(r, [[-1.0, 0.0, 5.0], [0.0, -2.0, -3.0]]))

    def test_pivot_is_fixed_point(self):
        pivot = np.array([3.0, -1.0, 2.0])
        pts = np.array([pivot, pivot + [1.0, 0.0, 0.0]])
        r = rotate_points_about_axis(pts, pivot=pivot, axis=[0, 0, 1], degrees=90)
        self.assertTrue(np.allclose(r[0], pivot))               # point on the axis is fixed
        self.assertTrue(np.allclose(r[1], pivot + [0.0, 1.0, 0.0]))  # +x -> +y about z

    def test_full_turn_is_identity(self):
        pts = np.random.default_rng(0).normal(size=(6, 3))
        r = rotate_points_about_axis(pts, pivot=[0.1, 0.2, 0.3], axis=[1, 1, 1], degrees=360)
        self.assertTrue(np.allclose(r, pts, atol=1e-9))

    def test_input_not_mutated(self):
        pts = np.array([[1.0, 2.0, 3.0]])
        orig = pts.copy()
        rotate_points_about_axis(pts, pivot=np.zeros(3), axis=[0, 1, 0], degrees=45)
        self.assertTrue(np.array_equal(pts, orig))


class TestPdbReplaceCoords(unittest.TestCase):
    def _write(self, text):
        fd, path = tempfile.mkstemp(suffix='.pdb')
        os.close(fd)
        with open(path, 'w') as fh:
            fh.write(text)
        return path

    def test_replaces_only_mapped_serials_and_preserves_columns(self):
        # two ATOM records; columns follow the fixed PDB format
        src = (
            "ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00      PROT N\n"
            "ATOM      2  CA  ALA A   1       4.000   5.000   6.000  1.00  0.00      PROT C\n"
            "END\n"
        )
        src_path = self._write(src)
        out_path = src_path + '.out.pdb'
        # coords array indexed by row; map serial 2 -> row 0 (serial 1 left untouched)
        coords = np.array([[7.123, 8.456, 9.789]])
        pdb_replace_coords(src_path, out_path, coords, {2: 0})
        with open(out_path) as fh:
            lines = fh.readlines()
        # serial 1 unchanged
        self.assertIn("1.000   2.000   3.000", lines[0])
        # serial 2 replaced, formatted %8.3f, other columns intact
        self.assertIn("   7.123   8.456   9.789", lines[1])
        self.assertTrue(lines[1].startswith("ATOM      2  CA  ALA A   1"))
        self.assertIn("PROT C", lines[1])
        self.assertEqual(lines[2].strip(), "END")
        for p in (src_path, out_path):
            os.remove(p)


if __name__ == '__main__':
    unittest.main()
