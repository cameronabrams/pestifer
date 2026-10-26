# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Unit tests for :mod:`pestifer.molecule.quilt`, the patch tiler.

Pure numpy on synthetic patches, so they run in CI.  A synthetic lipid is a few atoms: a
"phospholipid" is a vertical chain P -> C1 -> C2 reaching toward the midplane, and a "sterol" is
an O3 hydroxyl above a flat six-membered ring whose plane is horizontal, so a vertical bond can
thread it.

What is pinned, and why each one matters:

- Conflicts are found **across the periodic seam**.  A piercing whose two halves are a box length
  apart in raw coordinates is exactly the case ``RingChecker._check_fast`` once missed (see its
  docstring); a seam check that missed it would let every seam threading through.
- A box of whole tiles keeps every lipid and removes none -- the in-tile contacts are the patch's
  own, so a finding there would be a false alarm.
- After equalization the leaflet ratio is the patches' stress-free ratio to within a lipid.
- Leaflet membership comes from the headgroup reference atom, not the center of mass.
"""

import os
import tempfile
import unittest

import numpy as np

from pestifer.molecule.quilt import (Copy, Patch, PatchLipid, build_quilt, find_conflicts,
                                     resolve_conflicts, seams_of, tile_leaflet, write_quilt_pdb)

RING_R = 1.4


def _phospholipid(x, y, upper, tilt=0.0):
    """``tilt`` moves the tail laterally away from the phosphorus, as a real lipid's does, so its
    centroid and its headgroup are in different places."""
    s = 1.0 if upper else -1.0
    names = ['P', 'C1', 'C2']
    xyz = [(x, y, s * 18.0), (x + 0.5 * tilt, y, s * 10.0), (x + tilt, y, s * 3.0)]
    bonds = [(0, 1), (1, 2)]
    return 'POPC', names, xyz, bonds, 0


def _sterol(x, y, upper, z_ring=8.0, z_oh=14.0):
    s = 1.0 if upper else -1.0
    names = ['O3'] + [f'R{k}' for k in range(1, 7)]
    t = np.linspace(0, 2 * np.pi, 7)[:-1]
    xyz = [(x, y, s * z_oh)] + [(x + RING_R * np.cos(a), y + RING_R * np.sin(a), s * z_ring) for a in t]
    bonds = [(0, 1)] + [(k, k % 6 + 1) for k in range(1, 7)]
    return 'CHL1', names, xyz, bonds, 0


def _patch(molecules, cell, origin=None):
    """A :class:`Patch` from ``(resname, names, xyz, local_bonds, ref_local)`` tuples."""
    names, resnames, xyz, bonds, lipids = [], [], [], [], []
    for resname, nm, xs, bs, ref in molecules:
        base = len(names)
        rows = np.arange(base, base + len(nm))
        names += nm
        resnames += [resname] * len(nm)
        xyz += list(xs)
        bonds += [(base + a, base + b) for a, b in bs]
        lipids.append(PatchLipid(rows=rows, ref=base + ref, resname=resname))
    cell = np.asarray(cell, float)
    return Patch(names=np.array(names, dtype=object), resnames=np.array(resnames, dtype=object),
                 xyz=np.array(xyz, float), bonds=np.array(bonds, int).reshape(-1, 2),
                 lipids=lipids, cell=cell,
                 origin=(0.5 * cell if origin is None else np.asarray(origin, float)))


def _lattice_patch(n=4, spacing=8.0, sterol_every=0):
    """An n x n square lattice per leaflet, centered cell; every ``sterol_every``-th site a sterol."""
    mols = []
    k = 0
    for upper in (True, False):
        for i in range(n):
            for j in range(n):
                x, y = (i + 0.5) * spacing, (j + 0.5) * spacing
                if upper is False:
                    x, y = x + 0.25 * spacing, y + 0.25 * spacing   # leaflets need not register
                make = _sterol if sterol_every and k % sterol_every == 0 else _phospholipid
                mols.append(make(x, y, upper))
                k += 1
    L = n * spacing
    return _patch(mols, (L, L, 90.0))


def _random_patch(n=100, L=65.0, seed=1, tilt=0.0):
    """``n`` phospholipids per leaflet at random positions: the lateral granularity of a real
    relaxed patch, which a coarse lattice does not have (on an 8 A lattice the crop count jumps
    by whole rows as the box grows)."""
    rng = np.random.default_rng(seed)
    mols = [_phospholipid(*rng.uniform(0, L, 2), upper, tilt=tilt)
            for upper in (True, False) for _ in range(n)]
    return _patch(mols, (L, L, 90.0))


class TestRegister(unittest.TestCase):

    def test_wraps_whole_lipids_by_centroid_and_centers_the_midplane(self):
        P = _lattice_patch(n=2)
        P.xyz[:, 0] += 3 * P.cell[0] + 1.0     # unwrapped, as NAMD writes them
        P.xyz[:, 2] += 7.0
        P.register()
        for i in range(len(P.lipids)):
            c = P.centroid(i)
            self.assertTrue(0.0 <= c[0] < P.cell[0] and 0.0 <= c[1] < P.cell[1])
        refs = np.array([P.ref_point(i)[2] for i in range(len(P.lipids))])
        self.assertAlmostEqual(0.5 * (refs[refs > 0].mean() + refs[refs < 0].mean()), 0.0, places=6)

    def test_a_sunk_sterol_stays_in_the_leaflet_its_hydroxyl_points_out_of(self):
        """Its center of mass is below the midplane; its hydroxyl is above.  Membership by the
        hydroxyl is the rule; by mass it would be counted into the wrong leaflet."""
        mols = [_phospholipid(4, 4, True), _phospholipid(4, 4, False),
                _phospholipid(12, 12, True), _phospholipid(12, 12, False),
                _sterol(8, 8, True, z_ring=-6.0, z_oh=3.0)]
        P = _patch(mols, (16.0, 16.0, 90.0)).register()
        sterol = P.lipids[-1]
        self.assertLess(P.centroid(len(P.lipids) - 1)[2], 0.0)
        self.assertEqual(sterol.leaflet, 'upper')

    def test_rings_are_found_and_ordered_around_the_cycle(self):
        P = _patch([_sterol(5, 5, True), _phospholipid(9, 9, False)], (16.0, 16.0, 90.0)).register()
        self.assertEqual(len(P.rings[0]), 1)
        ring = P.xyz[P.rings[0][0]]
        steps = np.linalg.norm(np.roll(ring, -1, axis=0) - ring, axis=1)
        self.assertTrue(np.allclose(steps, RING_R, atol=1e-6), 'ring atoms are not in cycle order')
        self.assertEqual(P.rings[1], [])


class TestConflictsAcrossTheSeam(unittest.TestCase):
    """A sterol ring sitting just inside the far x edge and a phospholipid bond just inside the
    near one: a box length apart in raw coordinates, 0.2 A apart through the periodic seam."""

    L = 40.0

    def _copies(self, bond_x):
        P = _patch([_sterol(0.0, 0.0, True), _phospholipid(0.0, 0.0, True)], (8.0, 8.0, 90.0))
        P._index_topology()     # positions are set by hand below, so not register()
        ring = Copy(P, 0, np.array([self.L - 0.1, 20.0, 0.0]), 'upper')
        bond = Copy(P, 1, np.array([bond_x, 20.0, 0.0]), 'upper')
        return [ring, bond]

    def test_a_threading_through_the_periodic_seam_is_found(self):
        copies = self._copies(bond_x=0.1)
        kinds = [c.kind for c in find_conflicts(copies, self.L, self.L, clash_distance=0.3)]
        self.assertIn('pierce', kinds)

    def test_negative_control_the_same_pair_apart_is_not_a_threading(self):
        copies = self._copies(bond_x=0.5 * self.L)
        self.assertEqual(find_conflicts(copies, self.L, self.L, clash_distance=0.3), [])

    def test_a_clash_through_the_periodic_seam_is_found(self):
        P = _patch([_phospholipid(0, 0, True)], (8.0, 8.0, 90.0))
        P._index_topology()
        a = Copy(P, 0, np.array([self.L - 0.05, 10.0, 0.0]), 'upper')
        b = Copy(P, 0, np.array([0.05, 10.0, 0.0]), 'upper')
        far = Copy(P, 0, np.array([20.0, 10.0, 0.0]), 'upper')
        self.assertEqual([c.kind for c in find_conflicts([a, b], self.L, self.L, 0.3)][:1], ['clash'])
        self.assertEqual(find_conflicts([a, far], self.L, self.L, 0.3), [])

    def test_the_lipid_nearer_a_seam_is_the_one_removed(self):
        P = _patch([_phospholipid(0, 0, True)], (8.0, 8.0, 90.0))
        P._index_topology()
        near = Copy(P, 0, np.array([0.5, 20.0, 0.0]), 'upper')
        less_near = Copy(P, 0, np.array([self.L - 1.5, 20.0, 0.0]), 'upper')
        conflicts = find_conflicts([near, less_near], self.L, self.L, clash_distance=2.5)
        self.assertTrue(conflicts)
        gone = resolve_conflicts([near, less_near], conflicts, self.L, self.L, (True, True))
        self.assertEqual(gone, {0})


class TestTiling(unittest.TestCase):

    def test_whole_tiles_keep_every_lipid_and_have_no_seam(self):
        P = _lattice_patch(n=4, sterol_every=3).register()
        q = build_quilt(P, P, 3 * P.cell[0], 2 * P.cell[1], gap=2.0)
        self.assertEqual(q.seams, (False, False))
        self.assertEqual((q.n_upper, q.n_lower), (16 * 6, 16 * 6))
        self.assertEqual(q.removed, [])
        self.assertAlmostEqual(q.fill, 1.0)

    def test_whole_tiles_of_tilted_lipids_keep_every_lipid(self):
        """Tilted lipids put a headgroup and a centroid on opposite sides of a cell edge.  The
        wrap and the crop must use the same point, or such a lipid loses a copy at the edge of
        every whole-tile box."""
        P = _random_patch(tilt=4.0).register()
        q = build_quilt(P, P, 2 * P.cell[0], 2 * P.cell[1], gap=2.0)
        self.assertEqual((q.n_upper, q.n_lower), (400, 400))

    def test_a_seam_in_the_lower_leaflet_gaps_the_upper_too(self):
        """The upper patch tiles the box exactly; the lower does not.  If only the lower paid the
        gap, the two leaflets would lose different fractions of their area to it."""
        A = _random_patch(L=40.0, seed=2).register()
        B = _random_patch(L=37.0, seed=3).register()
        Lx = Ly = 3 * A.cell[0]
        gap = 4.0
        q = build_quilt(A, B, Lx, Ly, gap=gap)
        self.assertEqual(q.seams, (True, True))
        up = np.array([c.centroid for c in q.copies if c.leaflet == 'upper'])
        self.assertTrue((up[:, :2] >= gap / 2).all() and (up[:, :2] < Lx - gap / 2).all(),
                        'the upper leaflet was not gapped at the lower leaflet\'s seam')

    def test_a_seam_in_either_leaflet_gaps_both(self):
        """A gap paid by one leaflet only is a built-in error in the leaflet ratio."""
        A = _lattice_patch(n=4).register()
        B = _lattice_patch(n=4, spacing=7.0).register()
        self.assertEqual(seams_of([A, B], 2 * A.cell[0], 2 * A.cell[1]), (True, True))
        self.assertEqual(seams_of([A, A], 2 * A.cell[0], 2 * A.cell[1]), (False, False))

    def test_crop_keeps_lipids_by_centroid_inside_the_gap(self):
        P = _lattice_patch(n=4).register()
        Lx = Ly = 2.5 * P.cell[0]
        copies = tile_leaflet(P, 'upper', Lx, Ly, 4.0, (True, True))
        cen = np.array([c.centroid for c in copies])
        self.assertTrue(((cen[:, :2] >= 2.0) & (cen[:, :2] < Lx - 2.0)).all())

    def test_equalized_leaflets_carry_the_stress_free_ratio(self):
        """Two patches of different density: after cropping and equalizing, each leaflet holds
        the same fraction of the count its patch's density puts in the box, to within a lipid."""
        A = _lattice_patch(n=5, spacing=7.0).register()     # denser upper
        B = _lattice_patch(n=4, spacing=8.5).register()
        for Lx, Ly in ((100.0, 100.0), (83.0, 61.0), (120.0, 47.0)):
            with self.subTest(box=(Lx, Ly)):
                q = build_quilt(A, B, Lx, Ly, gap=2.0)
                fu = q.n_upper / q.stress_free['upper']
                fl = q.n_lower / q.stress_free['lower']
                self.assertLessEqual(abs(fu - fl) * min(q.stress_free.values()), 1.0 + 1e-9)

    def test_expand_relaxes_to_about_the_requested_area(self):
        P = _random_patch().register()
        Lx = Ly = 2.5 * P.cell[0]
        q = build_quilt(P, P, Lx, Ly, gap=4.0, expand=True)
        self.assertGreater(q.Lx, Lx)
        self.assertAlmostEqual(q.fill * q.Lx * q.Ly / (Lx * Ly), 1.0, delta=0.05)

    def test_box_height_gives_each_chamber_its_own_patch_half(self):
        A = _lattice_patch(n=4).register()
        B = _lattice_patch(n=4).register()
        B.cell[2] = 70.0
        q = build_quilt(A, B, 2 * A.cell[0], 2 * A.cell[1], gap=2.0)
        top = q.origin[2] + 0.5 * q.Lz
        bottom = q.origin[2] - 0.5 * q.Lz
        self.assertAlmostEqual(top, 0.5 * A.cell[2])
        self.assertAlmostEqual(bottom, -0.5 * B.cell[2])


class TestWriteQuiltPdb(unittest.TestCase):

    def test_layout_matches_the_grid_packer(self):
        mols = [_phospholipid(4, 4, True), _phospholipid(4, 4, False)]
        P = _patch(mols, (8.0, 8.0, 90.0)).register()
        P.lipids[0].resname = 'SAPI24'      # a six-character CHARMM name
        q = build_quilt(P, P, 16.0, 16.0, gap=2.0)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'q.pdb')
            write_quilt_pdb(q, path)
            lines = [ln for ln in open(path) if ln.startswith('ATOM')]
        self.assertEqual(len(lines), sum(len(c.patch.lipids[c.lipid].rows) for c in q.copies))
        resids = {int(ln[22:26]) for ln in lines}
        self.assertEqual(len(resids), len(q.copies), 'every lipid needs its own resid')
        xyz = np.array([(float(ln[30:38]), float(ln[38:46]), float(ln[46:54])) for ln in lines])
        want = np.concatenate([c.xyz for c in q.copies])
        self.assertTrue(np.allclose(xyz, want, atol=1e-3), 'coordinates are not in columns 31-54')
        self.assertEqual({ln[17:21] for ln in lines}, {'SAPI', 'POPC'})
        self.assertEqual({ln[21] for ln in lines}, {'A'})


    def test_a_long_name_is_written_as_its_placeholder(self):
        """Cut to four columns, SAPI24 reads as SAPI -- a different residue psfgen builds
        without complaint.  With an alias the writer puts the placeholder psfgen maps back."""
        mols = [_phospholipid(4, 4, True), _phospholipid(4, 4, False)]
        P = _patch(mols, (8.0, 8.0, 90.0)).register()
        P.lipids[0].resname = 'SAPI24'
        q = build_quilt(P, P, 16.0, 16.0, gap=2.0)
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, 'q.pdb')
            write_quilt_pdb(q, path, aliases={'SAPI24': 'Q001'})
            names = {ln[17:21] for ln in open(path) if ln.startswith('ATOM')}
        self.assertEqual(names, {'Q001', 'POPC'})


if __name__ == '__main__':
    unittest.main()
