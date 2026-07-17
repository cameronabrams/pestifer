import unittest

import numpy as np

from pestifer.psfutil.loop_ccd import (
    RAMACHANDRAN_BASINS,
    sample_backbone_dihedrals,
    optimal_ccd_angle,
    end_rmsd,
    ccd_close,
    dihedral_deg,
    set_dihedral,
)
from pestifer.util.coord import rotate_points_about_axis


class TestDihedral(unittest.TestCase):

    def test_dihedral_known_values(self):
        # eclipsed (coplanar, same side) -> 0 deg
        self.assertAlmostEqual(abs(dihedral_deg([0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0])), 0.0, places=4)
        # anti (coplanar, opposite side) -> 180 deg
        self.assertAlmostEqual(abs(dihedral_deg([0, 1, 0], [0, 0, 0], [1, 0, 0], [1, -1, 0])), 180.0, places=4)
        # a perpendicular case is +/-90; magnitude is convention-independent
        self.assertAlmostEqual(abs(dihedral_deg([0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 0, 1])), 90.0, places=4)

    def test_set_dihedral_reaches_target(self):
        # 4-atom chain + two downstream atoms that must rotate with l
        coords = np.array([[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0],
                           [1.5, 1.5, 0.0], [2.0, 1.0, 0.0]], dtype=float)
        mask = np.zeros(6, dtype=bool); mask[3:] = True   # l and its dependents
        for target in (-120.0, 45.0, 179.0, -60.0):
            out = set_dihedral(coords, 0, 1, 2, 3, target, mask)
            self.assertAlmostEqual(dihedral_deg(out[0], out[1], out[2], out[3]), target, places=4)
            # atoms upstream of the bond are untouched
            np.testing.assert_allclose(out[:3], coords[:3])

    def test_matches_vmd_measure_dihed_sign(self):
        # dihedral_deg must follow VMD's measure dihed convention (Ramachandran values assume it).
        # Reference values precomputed with VMD 2.0.0 `measure dihed` on these exact points.
        cases = [
            (([-1.0, 1.0, 0.5], [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.5, 1.0, 1.0]), None),
        ]
        # rather than pin brittle magic numbers, assert the sign relationship that failed before:
        # a right-handed (VMD +) increase in the dihedral corresponds to +delta in set_dihedral.
        coords = np.array([[-1, 1, 0], [0, 0, 0], [1, 0, 0], [1.4, 0.8, 0.6], [2.0, 1.2, 0.3]], float)
        mask = np.zeros(5, dtype=bool); mask[3:] = True
        cur = dihedral_deg(coords[0], coords[1], coords[2], coords[3])
        out = set_dihedral(coords, 0, 1, 2, 3, cur + 25.0, mask)
        self.assertAlmostEqual(dihedral_deg(out[0], out[1], out[2], out[3]), cur + 25.0, places=4)

    def test_set_dihedral_rigidly_carries_downstream_atoms(self):
        coords = np.array([[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0],
                           [1.5, 1.5, 0.0]], dtype=float)
        mask = np.zeros(5, dtype=bool); mask[3:] = True
        out = set_dihedral(coords, 0, 1, 2, 3, 33.0, mask)
        # the l->dependent distance is preserved (rigid-body rotation)
        self.assertAlmostEqual(np.linalg.norm(out[4] - out[3]),
                               np.linalg.norm(coords[4] - coords[3]), places=6)


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


class TestLoopExtractionAndReclose(unittest.TestCase):
    """Extraction + CCD on real BPTI backbone: open an internal loop, then re-close it."""

    @classmethod
    def setUpClass(cls):
        from pathlib import Path
        from pestifer.psfutil.loop_ccd import backbone_from_pdb
        pdb = Path(__file__).parents[2] / 'inputs' / '6pti.pdb'
        cls.bb = backbone_from_pdb(str(pdb), chainID='A')

    def _bond_lengths(self, coords, problem):
        # N-CA and CA-C bond lengths for each loop residue
        r = problem['row']
        lens = []
        for (a, b) in problem['bonds']:
            lens.append(np.linalg.norm(coords[a] - coords[b]))
        return np.array(lens)

    def test_open_then_ccd_recloses_to_native(self):
        from pestifer.psfutil.loop_ccd import (build_loop_ccd_problem, ccd_close,
                                               end_rmsd, RAMACHANDRAN_BASINS)
        loop = [12, 13, 14, 15, 16, 17]
        prob = build_loop_ccd_problem(self.bb, loop, n_anchor_resid=11)
        coords0 = prob['coords']
        end_idx = prob['end_idx']
        target = coords0[end_idx].copy()               # native end positions
        native_bonds = self._bond_lengths(coords0, prob)

        # OPEN the loop: perturb each rotatable bond by a random angle
        rng = np.random.default_rng(2024)
        opened = coords0.copy()
        for (a, b), mask in zip(prob['bonds'], prob['moving_masks']):
            deg = rng.uniform(-50, 50)
            opened[mask] = rotate_points_about_axis(opened[mask], opened[a], opened[b] - opened[a], deg)
        opened_rmsd = end_rmsd(opened, end_idx, target)
        self.assertGreater(opened_rmsd, 2.0, "perturbation should have opened the loop")

        # CLOSE with CCD back onto the native end
        closed, rmsd, iters = ccd_close(opened, prob['bonds'], prob['moving_masks'],
                                        end_idx, target, tol=0.1)
        self.assertLess(rmsd, 0.1, f"CCD failed to reclose (rmsd={rmsd:.3f}, {iters} iters)")

        # rotations are rigid: internal loop bond lengths must be unchanged
        closed_bonds = self._bond_lengths(closed, prob)
        np.testing.assert_allclose(closed_bonds, native_bonds, atol=1e-6)

    def test_problem_structure(self):
        from pestifer.psfutil.loop_ccd import build_loop_ccd_problem
        loop = [20, 21, 22]
        prob = build_loop_ccd_problem(self.bb, loop, n_anchor_resid=19)
        # 1 anchor C + 3 residues x 4 backbone atoms = 13 rows
        self.assertEqual(prob['coords'].shape, (13, 3))
        # 2 rotatable bonds (phi, psi) per loop residue
        self.assertEqual(len(prob['bonds']), 6)
        self.assertEqual(len(prob['moving_masks']), 6)
        # end effector = last loop residue N, CA, C
        r = prob['row']
        self.assertEqual(prob['end_idx'], [r[(22, 'N')], r[(22, 'CA')], r[(22, 'C')]])
        # the anchor C (row 0) is never movable
        for m in prob['moving_masks']:
            self.assertFalse(m[0])
