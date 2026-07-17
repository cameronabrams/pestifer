import unittest

import numpy as np

from pestifer.psfutil.loop_ccd import (
    RAMACHANDRAN_BASINS,
    sample_backbone_dihedrals,
    optimal_ccd_angle,
    end_rmsd,
    ccd_close,
)
from pestifer.util.coord import rotate_points_about_axis


class TestOptimalCCDAngle(unittest.TestCase):
    """theta* must be the closed-form minimizer of sum_k |R(theta) m_k - t_k|^2."""

    @staticmethod
    def _objective(moving, target, pivot, axis, deg):
        rot = rotate_points_about_axis(moving, pivot, axis, deg)
        d = rot - target
        return float(np.sum(d * d))

    def test_matches_bruteforce_minimum(self):
        rng = np.random.default_rng(12345)
        for _ in range(25):
            K = rng.integers(1, 5)
            moving = rng.normal(size=(K, 3))
            target = rng.normal(size=(K, 3))
            pivot = rng.normal(size=3)
            axis = rng.normal(size=3)
            theta = optimal_ccd_angle(moving, target, pivot, axis)
            # brute-force scan
            grid = np.linspace(-180.0, 180.0, 3601)  # 0.1 deg resolution
            vals = [self._objective(moving, target, pivot, axis, g) for g in grid]
            best = grid[int(np.argmin(vals))]
            # angles are circular; compare objective values instead of raw angle
            self.assertLessEqual(
                self._objective(moving, target, pivot, axis, theta),
                self._objective(moving, target, pivot, axis, best) + 1e-6,
            )

    def test_zero_when_moving_equals_target(self):
        # already on target -> no rotation improves it
        pts = np.array([[1.0, 2.0, 0.5], [0.0, -1.0, 2.0]])
        theta = optimal_ccd_angle(pts, pts.copy(), pivot=[0, 0, 0], axis=[0, 0, 1])
        # objective is already minimal; rotating by theta must not increase it
        self.assertAlmostEqual(
            self._objective(pts, pts, [0, 0, 0], [0, 0, 1], theta), 0.0, places=6)

    def test_recovers_a_known_rotation(self):
        rng = np.random.default_rng(7)
        moving = rng.normal(size=(4, 3))
        pivot = np.array([0.3, -0.2, 0.1])
        axis = np.array([0.0, 0.0, 1.0])
        applied = 37.0
        target = rotate_points_about_axis(moving, pivot, axis, applied)
        theta = optimal_ccd_angle(moving, target, pivot, axis)
        self.assertAlmostEqual(theta, applied, places=4)

    def test_zero_axis_returns_zero(self):
        self.assertEqual(optimal_ccd_angle([[1, 0, 0]], [[0, 1, 0]], [0, 0, 0], [0, 0, 0]), 0.0)


class TestCCDClose(unittest.TestCase):
    """ccd_close should drive the end effector onto a reachable target."""

    @staticmethod
    def _chain(M):
        # a gently bent chain in 3D so every bond axis moves downstream atoms
        t = np.linspace(0, 1, M)
        return np.stack([t * 10.0, np.sin(t * 3.0) * 2.0, np.cos(t * 2.0)], axis=1)

    @staticmethod
    def _bonds_and_masks(M):
        # bond i is (i, i+1); rotating it moves every atom strictly downstream of atom i+1
        bonds, masks = [], []
        for i in range(M - 1):
            bonds.append((i, i + 1))
            m = np.zeros(M, dtype=bool)
            m[i + 2:] = True
            masks.append(m)
        return bonds, masks

    def test_converges_to_reachable_target(self):
        M = 10
        coords = self._chain(M)
        bonds, masks = self._bonds_and_masks(M)
        end_idx = [M - 1]
        # build a provably reachable target by applying random rotations to the chain
        rng = np.random.default_rng(2024)
        moved = coords.copy()
        for (a, b), m in zip(bonds, masks):
            deg = rng.uniform(-40, 40)
            moved[m] = rotate_points_about_axis(moved[m], moved[a], moved[b] - moved[a], deg)
        target = moved[end_idx]
        closed, rmsd, iters = ccd_close(coords, bonds, masks, end_idx, target, tol=0.05)
        self.assertLess(rmsd, 0.05, f"CCD did not close (rmsd={rmsd:.3f} after {iters} iters)")

    def test_does_not_move_atoms_upstream_of_all_bonds(self):
        # the first two atoms (0,1) are on/at the first bond and never in any moving mask
        M = 8
        coords = self._chain(M)
        bonds, masks = self._bonds_and_masks(M)
        target = coords[[M - 1]] + np.array([1.0, 1.0, 1.0])
        closed, _, _ = ccd_close(coords, bonds, masks, [M - 1], target, max_iters=50)
        np.testing.assert_allclose(closed[0], coords[0])
        np.testing.assert_allclose(closed[1], coords[1])

    def test_returns_copy_not_mutating_input(self):
        M = 6
        coords = self._chain(M)
        original = coords.copy()
        bonds, masks = self._bonds_and_masks(M)
        ccd_close(coords, bonds, masks, [M - 1], coords[[M - 1]] + 1.0, max_iters=20)
        np.testing.assert_array_equal(coords, original)

    def test_unreachable_target_terminates_without_error(self):
        # target far beyond the chain's reach: must stop (stall/max_iters), not hang or throw
        M = 6
        coords = self._chain(M)
        bonds, masks = self._bonds_and_masks(M)
        far = coords[[M - 1]] + np.array([100.0, 100.0, 100.0])
        closed, rmsd, iters = ccd_close(coords, bonds, masks, [M - 1], far, max_iters=100)
        self.assertGreater(rmsd, 1.0)   # obviously can't reach
        self.assertLessEqual(iters, 100)


class TestSampleBackboneDihedrals(unittest.TestCase):

    def test_deterministic_for_fixed_seed(self):
        a = sample_backbone_dihedrals(np.random.default_rng(42), 5)
        b = sample_backbone_dihedrals(np.random.default_rng(42), 5)
        np.testing.assert_array_equal(a, b)

    def test_shape_and_within_basin_neighborhood(self):
        out = sample_backbone_dihedrals(np.random.default_rng(1), 200, jitter_deg=8.0)
        self.assertEqual(out.shape, (200, 2))
        centers = np.array([[v[0], v[1]] for v in RAMACHANDRAN_BASINS.values()])
        # every sample must be near (< 6 sigma) at least one basin center
        for phi, psi in out:
            dmin = np.min(np.hypot(centers[:, 0] - phi, centers[:, 1] - psi))
            self.assertLess(dmin, 6 * 8.0 * np.sqrt(2))

    def test_weights_bias_toward_dominant_basins(self):
        # alphaR+beta together (0.75 weight) should dominate a large sample
        out = sample_backbone_dihedrals(np.random.default_rng(9), 4000)
        near_alphaR = np.sum(np.hypot(out[:, 0] + 63, out[:, 1] + 43) < 30)
        near_alphaL = np.sum(np.hypot(out[:, 0] - 60, out[:, 1] - 45) < 30)
        self.assertGreater(near_alphaR, near_alphaL * 3)   # 0.42 vs 0.05


if __name__ == '__main__':
    unittest.main()
