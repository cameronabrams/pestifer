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
