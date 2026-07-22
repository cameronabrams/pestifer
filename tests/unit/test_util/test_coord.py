import os
import tempfile
import numpy as np
import unittest

from pestifer.util.coord import (
    rotate_points_about_axis,
    pdb_replace_coords,
    standardize_pdb_columns,
    apply_tmat,
    kabsch,
    build_tmat,
    orient_peptide_fusion,
)


def _unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


class TestOrientPeptideFusion(unittest.TestCase):
    def _case(self, seed=0):
        rng = np.random.default_rng(seed)
        base_C = rng.normal(size=3)
        base_CA = base_C + _unit(rng.normal(size=3)) * 1.5
        base_O = base_C + _unit(rng.normal(size=3)) * 1.23
        donor_N = rng.normal(size=3) + 25
        donor_CA = donor_N + _unit(rng.normal(size=3)) * 1.47
        donor = np.vstack([donor_N, donor_CA, donor_N + rng.normal(size=(5, 3))])
        return base_C, base_CA, base_O, donor_N, donor_CA, donor

    def test_N_lands_at_ideal_amide_site(self):
        bC, bCA, bO, dN, dCA, donor = self._case(1)
        out = orient_peptide_fusion(donor, dN, dCA, bC, bCA, bO)
        # donor N (row 0) is 1.33 A from the base carbonyl C
        self.assertAlmostEqual(np.linalg.norm(out[0] - bC), 1.33, places=6)

    def test_first_NCA_aligns_with_CN_bond(self):
        bC, bCA, bO, dN, dCA, donor = self._case(2)
        out = orient_peptide_fusion(donor, dN, dCA, bC, bCA, bO)
        n_ca = _unit(out[1] - out[0])
        c_n = _unit(out[0] - bC)
        self.assertTrue(np.allclose(n_ca, c_n, atol=1e-9))

    def test_matches_tcl_reference_math(self):
        bC, bCA, bO, dN, dCA, donor = self._case(3)
        # independent reimplementation of the VMD/Tcl vecmath
        dirN = _unit(-1.0 * (_unit(bCA - bC) + _unit(bO - bC)))
        Nideal = bC + 1.33 * dirN
        moved = donor + (Nideal - dN)
        u = _unit((dCA + (Nideal - dN)) - Nideal)
        v = _unit(Nideal - bC)
        axis = np.cross(u, v)
        c = np.clip(np.dot(u, v), -1, 1)
        ref = rotate_points_about_axis(moved, Nideal, axis, np.degrees(np.arccos(c)))
        out = orient_peptide_fusion(donor, dN, dCA, bC, bCA, bO)
        self.assertTrue(np.allclose(out, ref, atol=1e-12))

    def test_input_not_mutated(self):
        bC, bCA, bO, dN, dCA, donor = self._case(4)
        orig = donor.copy()
        orient_peptide_fusion(donor, dN, dCA, bC, bCA, bO)
        self.assertTrue(np.array_equal(donor, orig))

    def test_collinear_NCA_and_CN_is_translation_only(self):
        # donor N->CA already parallel to the C->N bond: axis ~ 0, no rotation
        bC = np.zeros(3)
        bCA = np.array([1.5, 0.0, 0.0])
        bO = np.array([-0.5, 1.1, 0.0])
        dirN = _unit(-1.0 * (_unit(bCA - bC) + _unit(bO - bC)))
        Nideal = bC + 1.33 * dirN
        dN = np.array([10.0, 10.0, 10.0])
        dCA = dN + _unit(Nideal - bC) * 1.47      # N->CA parallel to C->N
        donor = np.vstack([dN, dCA])
        out = orient_peptide_fusion(donor, dN, dCA, bC, bCA, bO)
        self.assertTrue(np.allclose(out, donor + (Nideal - dN), atol=1e-9))


class TestApplyTmat(unittest.TestCase):
    def test_identity_leaves_points(self):
        pts = np.array([[1.0, 2.0, 3.0], [-4.0, 5.0, -6.0]])
        self.assertTrue(np.allclose(apply_tmat(np.identity(4), pts), pts))

    def test_pure_translation(self):
        tmat = build_tmat(np.identity(3), np.array([1.0, 2.0, 3.0]))
        pts = np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
        self.assertTrue(np.allclose(apply_tmat(tmat, pts), [[1.0, 2.0, 3.0], [2.0, 3.0, 4.0]]))

    def test_single_point_shape_preserved(self):
        tmat = build_tmat(np.identity(3), np.array([5.0, 0.0, 0.0]))
        out = apply_tmat(tmat, np.array([1.0, 1.0, 1.0]))
        self.assertEqual(out.shape, (3,))
        self.assertTrue(np.allclose(out, [6.0, 1.0, 1.0]))

    def test_rotation_then_translation_order(self):
        # 90 deg about z (x->y, y->-x) then translate by +x
        R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        tmat = build_tmat(R, np.array([10.0, 0.0, 0.0]))
        out = apply_tmat(tmat, np.array([[1.0, 0.0, 0.0]]))
        # R @ (1,0,0) = (0,1,0); + (10,0,0) = (10,1,0)
        self.assertTrue(np.allclose(out, [[10.0, 1.0, 0.0]]))

    def test_input_not_mutated(self):
        pts = np.array([[1.0, 2.0, 3.0]])
        orig = pts.copy()
        apply_tmat(build_tmat(np.identity(3), np.array([1.0, 1.0, 1.0])), pts)
        self.assertTrue(np.array_equal(pts, orig))


class TestKabsch(unittest.TestCase):
    def _random(self, n=12, seed=0):
        return np.random.default_rng(seed).normal(size=(n, 3))

    def test_recovers_known_rigid_transform(self):
        P = self._random()
        R = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
        t = np.array([3.0, -2.0, 1.0])
        Q = P @ R.T + t
        Rf, tf, rmsd = kabsch(P, Q)
        self.assertTrue(np.allclose(Rf, R, atol=1e-9))
        self.assertTrue(np.allclose(tf, t, atol=1e-9))
        self.assertAlmostEqual(rmsd, 0.0, places=9)

    def test_result_is_proper_rotation(self):
        Rf, _, _ = kabsch(self._random(seed=1), self._random(seed=2))
        self.assertAlmostEqual(np.linalg.det(Rf), 1.0, places=9)
        self.assertTrue(np.allclose(Rf @ Rf.T, np.identity(3), atol=1e-9))

    def test_roundtrip_via_apply_tmat(self):
        P = self._random(seed=3)
        R = rotate_points_about_axis(np.identity(3), np.zeros(3), [1, 1, 1], 37.0)
        Q = P @ R.T + np.array([1.0, 2.0, 3.0])
        Rf, tf, _ = kabsch(P, Q)
        moved = apply_tmat(build_tmat(Rf, tf), P)
        self.assertTrue(np.allclose(moved, Q, atol=1e-9))

    def test_reflection_not_introduced_for_coplanar_points(self):
        # points in the z=0 plane; a naive SVD fit can flip handedness
        P = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])
        Q = P.copy()
        Rf, _, _ = kabsch(P, Q)
        self.assertAlmostEqual(np.linalg.det(Rf), 1.0, places=9)

    def test_rejects_mismatched_shapes(self):
        with self.assertRaises(ValueError):
            kabsch(self._random(4), self._random(5))


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


class TestStandardizePdbColumns(unittest.TestCase):
    def _write(self, text):
        fd, path = tempfile.mkstemp(suffix='.pdb')
        with os.fdopen(fd, 'w') as fh:
            fh.write(text)
        return path

    def test_reanchors_wide_resname_coords_and_leaves_standard_lines(self):
        # first line: standard protein record (coords already at cols 31-54)
        # second line: psfgen wide record -- 6-char carb resname BGLCNA shifts coords right by ~3
        src = (
            "ATOM      1  N   ALA A 512     193.912 186.468 172.509  0.00  0.00      A     \n"
            "ATOM  30841  C4  BGLCNAA 665     218.343 223.129 159.086  1.00  0.00      AG01 C\n"
            "END\n"
        )
        path = self._write(src)
        n = standardize_pdb_columns(path)
        self.assertEqual(n, 1)  # only the wide glycan line moved
        with open(path) as fh:
            lines = fh.readlines()
        # protein line unchanged
        self.assertEqual(lines[0].rstrip(),
                         "ATOM      1  N   ALA A 512     193.912 186.468 172.509  0.00  0.00      A")
        # glycan coords now parse from the standard columns
        self.assertAlmostEqual(float(lines[1][30:38]), 218.343, places=3)
        self.assertAlmostEqual(float(lines[1][38:46]), 223.129, places=3)
        self.assertAlmostEqual(float(lines[1][46:54]), 159.086, places=3)
        # trailing occupancy/beta/segname preserved
        self.assertIn("1.00  0.00      AG01", lines[1])
        # idempotent: a second pass moves nothing
        self.assertEqual(standardize_pdb_columns(path), 0)
        os.remove(path)


if __name__ == '__main__':
    unittest.main()
