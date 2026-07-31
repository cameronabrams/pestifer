# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for the athermal (repulsive-only) MC conformer sampler engine.

These exercise the force-field-agnostic core on a synthetic all-trans alkane-like chain, so
they run fast and deterministically with no CHARMM/psfgen dependency.
"""
import networkx as nx
import numpy as np
import pytest

from pestifer.charmmff.athermal_mc import (
    MoleculeMC, RotatableBond, run_mc, radial_distances, build_exclusions, moving_set,
    build_lipid_mc, cylinder_radius_for_apl,
)


def _trans_chain(n=16, bond=1.53, angle_deg=112.0):
    """An all-trans (planar zig-zag) carbon backbone: n atoms, uniform bond length/angle.

    Returns coords (n,3), a linear-chain graph, and the per-atom radii.  The chain lies in the
    x-z plane extended along +z, mimicking the extended-rod starting conformer the real sampler
    must melt.
    """
    half = np.deg2rad(angle_deg) / 2.0
    dz = bond * np.sin(half)
    dx = bond * np.cos(half)
    coords = np.zeros((n, 3))
    for i in range(1, n):
        coords[i, 2] = coords[i - 1, 2] + dz
        coords[i, 0] = coords[i - 1, 0] + (dx if i % 2 else -dx)
    g = nx.path_graph(n)
    radii = np.full(n, 1.7 * 2.0 ** (-1.0 / 6.0))   # carbon-ish sigma/2
    return coords, g, radii


def _chain_mol(n=16, cylinder_radius=float('inf')):
    coords, g, radii = _trans_chain(n)
    excl = build_exclusions(g, order=2)
    # every backbone bond except the first (anchor) is a pivot; tip side is the higher indices
    rot = [RotatableBond(a=i, b=i + 1, moving=moving_set(g, i, i + 1)) for i in range(1, n - 1)]
    mol = MoleculeMC(coords=coords, radii=radii, rotatable=rot, exclusions=excl,
                     axis_point=coords[0], axis_dir=np.array([0.0, 0.0, 1.0]),
                     cylinder_radius=cylinder_radius)
    return mol, g


def _bond_lengths(coords, g):
    return np.array([np.linalg.norm(coords[i] - coords[j]) for i, j in g.edges])


class TestAthermalMCEngine:

    def test_moving_set_linear_chain(self):
        _, g = _chain_mol(n=8)
        assert list(moving_set(g, 3, 4)) == [4, 5, 6, 7]
        assert list(moving_set(g, 0, 1)) == [1, 2, 3, 4, 5, 6, 7]

    def test_moving_set_rejects_ring_bond(self):
        g = nx.cycle_graph(6)
        with pytest.raises(ValueError):
            moving_set(g, 0, 1)

    def test_exclusions_cover_12_and_13(self):
        g = nx.path_graph(5)
        excl = build_exclusions(g, order=2)
        assert frozenset((0, 1)) in excl   # 1-2
        assert frozenset((0, 2)) in excl   # 1-3
        assert frozenset((0, 3)) not in excl   # 1-4 is allowed to clash

    def test_pivot_preserves_bonds_and_angles(self):
        """A dihedral pivot is a rigid rotation of the tip subtree, so every bond length and bond
        angle is invariant -- the whole point of torsion-only MC on fixed CHARMM geometry."""
        mol, g = _chain_mol(n=16)
        b0 = _bond_lengths(mol.coords, g)

        def angles(coords):
            out = []
            for k in range(1, len(coords) - 1):
                u = coords[k - 1] - coords[k]
                v = coords[k + 1] - coords[k]
                out.append(np.degrees(np.arccos(
                    np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))))
            return np.array(out)

        a0 = angles(mol.coords)
        samples = run_mc(mol, nsamples=5, n_equil=500, n_decorr=100, seed=7)
        for s in samples:
            assert np.allclose(_bond_lengths(s, g), b0, atol=1e-8)
            assert np.allclose(angles(s), a0, atol=1e-6)

    def test_no_overlaps_in_samples(self):
        """Every accepted conformer must be free of hard-sphere overlaps among non-excluded
        pairs -- the excluded-volume half of the athermal criterion."""
        mol, g = _chain_mol(n=16)
        samples = run_mc(mol, nsamples=6, n_equil=800, n_decorr=150, seed=3)
        for s in samples:
            for i in range(len(s)):
                for j in range(i + 1, len(s)):
                    if frozenset((i, j)) in mol.exclusions:
                        continue
                    d = np.linalg.norm(s[i] - s[j])
                    assert d >= mol.radii[i] + mol.radii[j] - 1e-9, \
                        f'overlap between {i},{j}: {d:.3f}'

    def test_confinement_keeps_atoms_in_cylinder(self):
        """With a tight cylinder, no confined atom in any sample may exceed the radius."""
        R = 4.0
        mol, g = _chain_mol(n=20, cylinder_radius=R)
        samples = run_mc(mol, nsamples=6, n_equil=1500, n_decorr=200, seed=11)
        for s in samples:
            rd = radial_distances(s[mol.confined], mol.axis_point, mol.axis_dir)
            assert rd.max() <= R + 1e-9

    def test_sampler_melts_the_rod(self):
        """Confined athermal MC must actually compress the extended rod: the axial span of the
        ensemble should fall well below the all-trans starting length (rods pack to low APL; the
        melt is what seeds fluid-like footprints)."""
        mol, g = _chain_mol(n=20, cylinder_radius=4.5)
        start_span = np.ptp(mol.coords[:, 2])
        samples = run_mc(mol, nsamples=8, n_equil=2500, n_decorr=250, seed=5)
        spans = np.array([np.ptp(s[:, 2]) for s in samples])
        assert spans.mean() < 0.9 * start_span

    def test_determinism(self):
        """Same seed -> identical trajectory (reproducible cached ensembles)."""
        mol1, _ = _chain_mol(n=16, cylinder_radius=5.0)
        mol2, _ = _chain_mol(n=16, cylinder_radius=5.0)
        s1 = run_mc(mol1, nsamples=4, n_equil=600, n_decorr=120, seed=42)
        s2 = run_mc(mol2, nsamples=4, n_equil=600, n_decorr=120, seed=42)
        for a, b in zip(s1, s2):
            assert np.array_equal(a, b)

    def test_no_rotatable_bonds_returns_input(self):
        coords, g, radii = _trans_chain(n=6)
        mol = MoleculeMC(coords=coords, radii=radii, rotatable=[])
        samples = run_mc(mol, nsamples=3, seed=1)
        assert len(samples) == 3
        for s in samples:
            assert np.array_equal(s, coords)


def _synthetic_diacyl():
    """A minimal two-tail 'lipid': P-O-glycerol, two ester carbonyls (each C=O), and two
    all-carbon tails.  The ester oxygens are what must stop the all-carbon tail-peeling, so the
    ester/glycerol/head bonds stay rigid and only the tail C-C bonds become rotatable."""
    #   idx: 0=P 1=O 2=C(gly) 3=C(=O)A 4=C 5=C 6=C(term) 7=C(=O)B 8=C 9=C 10=C(term)
    #        11=O(=A) 12=O(=B)
    elements = ['P', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'O', 'O']
    masses = [30.974, 15.999, 12.011, 12.011, 12.011, 12.011, 12.011,
              12.011, 12.011, 12.011, 12.011, 15.999, 15.999]
    bonds = [(0, 1), (1, 2), (2, 3), (3, 11), (3, 4), (4, 5), (5, 6),
             (2, 7), (7, 12), (7, 8), (8, 9), (9, 10)]
    coords = np.zeros((13, 3))
    coords[0] = (0.0, 0.0, 6.0)     # head P, highest z
    coords[1] = (0.0, 0.0, 5.0)
    coords[2] = (0.0, 0.0, 4.0)     # glycerol
    coords[3] = (-1.0, 0.0, 3.0); coords[11] = (-2.0, 0.0, 3.0)   # carbonyl A + =O
    coords[4] = (-1.0, 0.0, 2.0); coords[5] = (-1.0, 0.0, 1.0); coords[6] = (-1.0, 0.0, 0.0)
    coords[7] = (1.0, 0.0, 3.0);  coords[12] = (2.0, 0.0, 3.0)    # carbonyl B + =O
    coords[8] = (1.0, 0.0, 2.0);  coords[9] = (1.0, 0.0, 1.0); coords[10] = (1.0, 0.0, 0.0)
    rmin_half = np.full(13, 1.0)
    return coords, elements, masses, bonds, rmin_half


class TestBuildLipidMC:

    def test_selects_only_tail_torsions(self):
        coords, elements, masses, bonds, rmin_half = _synthetic_diacyl()
        mol = build_lipid_mc(coords, elements, masses, bonds,
                             head_indices=[0], tail_indices=[6, 10],
                             rmin_half=rmin_half, cylinder_radius=5.0, min_tip_heavy=2)
        pivots = {(rb.a, rb.b) for rb in mol.rotatable}
        # C(=O)-C and C-C tail torsions, but NOT terminal methyls, NOT ester/glycerol/head bonds
        assert pivots == {(3, 4), (4, 5), (7, 8), (8, 9)}

    def test_moving_sets_are_tip_subtrees(self):
        coords, elements, masses, bonds, rmin_half = _synthetic_diacyl()
        mol = build_lipid_mc(coords, elements, masses, bonds,
                             head_indices=[0], tail_indices=[6, 10],
                             rmin_half=rmin_half, cylinder_radius=5.0)
        mv = {(rb.a, rb.b): tuple(rb.moving) for rb in mol.rotatable}
        assert mv[(3, 4)] == (4, 5, 6)
        assert mv[(4, 5)] == (5, 6)
        assert mv[(7, 8)] == (8, 9, 10)

    def test_confines_tail_atoms_only(self):
        coords, elements, masses, bonds, rmin_half = _synthetic_diacyl()
        mol = build_lipid_mc(coords, elements, masses, bonds,
                             head_indices=[0], tail_indices=[6, 10],
                             rmin_half=rmin_half, cylinder_radius=5.0)
        confined = set(np.nonzero(mol.confined)[0].tolist())
        assert confined == {4, 5, 6, 8, 9, 10}   # both acyl tails, nothing in head/ester

    def test_radii_and_axis(self):
        coords, elements, masses, bonds, rmin_half = _synthetic_diacyl()
        mol = build_lipid_mc(coords, elements, masses, bonds,
                             head_indices=[0], tail_indices=[6, 10],
                             rmin_half=rmin_half, cylinder_radius=5.0, radius_scale=0.9)
        assert np.allclose(mol.radii, 0.9)
        # confinement axis is the membrane normal (+z by default), anchored on the tail bundle
        assert np.allclose(mol.axis_dir, [0.0, 0.0, 1.0])
        tail_xy = coords[[4, 5, 6, 8, 9, 10]].mean(axis=0)
        assert np.allclose(mol.axis_point, tail_xy)

    def test_elements_inferred_from_mass(self):
        coords, elements, masses, bonds, rmin_half = _synthetic_diacyl()
        mol = build_lipid_mc(coords, None, masses, bonds,
                             head_indices=[0], tail_indices=[6, 10],
                             rmin_half=rmin_half, cylinder_radius=5.0)
        assert {(rb.a, rb.b) for rb in mol.rotatable} == {(3, 4), (4, 5), (7, 8), (8, 9)}

    def test_run_mc_on_synthetic_lipid_preserves_geometry(self):
        coords, elements, masses, bonds, rmin_half = _synthetic_diacyl()
        mol = build_lipid_mc(coords, elements, masses, bonds,
                             head_indices=[0], tail_indices=[6, 10],
                             rmin_half=np.full(13, 0.6), cylinder_radius=6.0)
        b0 = np.array([np.linalg.norm(coords[i] - coords[j]) for i, j in bonds])
        samples = run_mc(mol, nsamples=4, n_equil=400, n_decorr=100, seed=9)
        for s in samples:
            bl = np.array([np.linalg.norm(s[i] - s[j]) for i, j in bonds])
            assert np.allclose(bl, b0, atol=1e-8)

    def test_cylinder_radius_for_apl(self):
        assert cylinder_radius_for_apl(60.0) == pytest.approx((60.0 / np.pi) ** 0.5, rel=1e-9)
        # ~8.7 A diameter for APL 60, per the design doc
        assert 2 * cylinder_radius_for_apl(60.0) == pytest.approx(8.74, abs=0.05)
