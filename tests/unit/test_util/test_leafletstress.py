"""Unit tests for pestifer.util.dcd and pestifer.util.leafletstress (no NAMD).

The geometry here decides which leaflet a slab's tension is charged to, and an error of one slab
moved example 17's Dgamma by 8-15 mN/m, so each rule is pinned against a case whose answer is
known by construction: a synthetic DCD with a bilayer of known midplane, and profiles whose
leaflet integrals can be done by hand.
"""
import struct

import numpy as np
import pytest

from pestifer.util.dcd import iter_frames, read_header
from pestifer.util.leafletstress import (
    bilayer_midplanes, leaflet_tensions, read_xst, slab_tensions, _periodic_cumulative)


def write_dcd(path, frames, cells, istart=100, nsavc=100):
    """A NAMD-shaped DCD: ``frames`` is a list of (natoms, 3) arrays, ``cells`` of (a, b, c)."""
    natoms = len(frames[0])

    def rec(b):
        return struct.pack('<i', len(b)) + b + struct.pack('<i', len(b))

    icntrl = [0] * 20
    icntrl[0], icntrl[1], icntrl[2] = len(frames), istart, nsavc
    icntrl[10], icntrl[19] = 1, 24                 # unit cell present; CHARMM version
    out = rec(b'CORD' + struct.pack('<20i', *icntrl))
    out += rec(struct.pack('<i', 1) + b'test'.ljust(80))
    out += rec(struct.pack('<i', natoms))
    for xyz, (a, b, c) in zip(frames, cells):
        out += rec(struct.pack('<6d', a, 90.0, b, 90.0, 90.0, c))
        for d in range(3):
            out += rec(np.asarray(xyz[:, d], dtype='<f4').tobytes())
    with open(path, 'wb') as f:
        f.write(out)


class TestDCDReader:
    def test_header_steps_cell_and_coordinates(self, tmp_path):
        rng = np.random.default_rng(0)
        frames = [rng.normal(size=(5, 3)) for _ in range(3)]
        p = tmp_path / 't.dcd'
        write_dcd(p, frames, [(10.0, 11.0, 12.0 + k) for k in range(3)], istart=600, nsavc=500)
        h = read_header(str(p))
        assert (h.natoms, h.nframes, h.has_cell) == (5, 3, True)
        # frame k is written at ISTART + k*NSAVC: firsttimestep + dcdfreq is the FIRST frame
        assert list(h.steps()) == [600, 1100, 1600]
        got = list(iter_frames(str(p), atoms=np.array([1, 3])))
        assert got[2][1] == (10.0, 11.0, 14.0)
        assert np.allclose(got[1][2], frames[1][[1, 3]], atol=1e-6)

    def test_the_frame_count_comes_from_the_file_not_nset(self, tmp_path):
        # a killed run leaves NSET stale; a truncated last frame must not be read
        p = tmp_path / 't.dcd'
        write_dcd(p, [np.zeros((4, 3))] * 3, [(10.0, 10.0, 10.0)] * 3)
        data = p.read_bytes()
        p.write_bytes(data[:-10])
        assert read_header(str(p)).nframes == 2


class TestPeriodicIntegration:
    def test_cumulative_wraps(self):
        g = np.array([1.0, 2.0, 3.0, 4.0])
        assert _periodic_cumulative(g, 2.5) == pytest.approx(1 + 2 + 1.5)
        assert _periodic_cumulative(g, 6.0) == pytest.approx(10 + 3)

    def test_leaflets_are_half_cells_on_either_side_of_the_midplane(self):
        # 8 slabs of 1 A in a cell from 0 to 8 (origin 4); midplane at z=5 -> slab coordinate 5
        g = np.arange(8, dtype=float)[None, :]
        lower, upper = leaflet_tensions(g, np.array([5.0]), np.array([4.0]), np.array([8.0]))
        assert lower[0] == pytest.approx(1 + 2 + 3 + 4)          # slabs 1..4
        assert upper[0] == pytest.approx(5 + 6 + 7 + 0)          # slabs 5,6,7 then wrap to 0
        assert lower[0] + upper[0] == pytest.approx(g.sum())

    def test_a_drifted_bilayer_gives_the_same_leaflet_tensions(self):
        """The reason for measuring the midplane: the same membrane shifted by 3 slabs."""
        base = np.array([0, 1, -5, 9, 0, 0, 7, -2, 4, 0], dtype=float)  # bilayer around slab 5
        lo0, up0 = leaflet_tensions(base[None], np.array([5.0]), np.array([5.0]),
                                    np.array([10.0]))
        shifted = np.roll(base, 3)
        lo1, up1 = leaflet_tensions(shifted[None], np.array([8.0]), np.array([5.0]),
                                    np.array([10.0]))
        assert (lo1[0], up1[0]) == pytest.approx((lo0[0], up0[0]))
        # and splitting the shifted one at the central slab gets it wrong
        assert shifted[5:].sum() - shifted[:5].sum() != pytest.approx(up0[0] - lo0[0])

    def test_a_fractional_midplane_splits_its_slab(self):
        g = np.ones((1, 4))
        lower, upper = leaflet_tensions(g, np.array([2.25]), np.array([2.0]), np.array([4.0]))
        assert lower[0] == pytest.approx(2.0) and upper[0] == pytest.approx(2.0)

    def test_slab_tension_units(self):
        # 100 bar of Pzz - PL over a 2 A slab = 0.01 * 100 * 2 = 2 mN/m
        p = np.zeros((1, 3, 3))
        p[0, :, 2] = 100.0
        assert np.allclose(slab_tensions(p, np.array([6.0])), 2.0)


class TestMidplane:
    def _bilayer(self, tmp_path, center, c=60.0, nframes=3, wrap=True):
        """Two leaflets of 'lipids': heads at center +/- 18, tails reaching center +/- 2."""
        rng = np.random.default_rng(1)
        n = 200
        z_up = center + rng.uniform(2.0, 18.0, n)
        z_lo = center - rng.uniform(2.0, 18.0, n)
        z = np.r_[z_up, z_lo]
        if wrap:                                   # wrap into a cell from 0 to c
            z = np.mod(z, c)
        xyz = np.zeros((len(z), 3))
        xyz[:, 2] = z
        p = tmp_path / 'b.dcd'
        write_dcd(p, [xyz] * nframes, [(40.0, 40.0, c)] * nframes)
        return str(p), np.arange(len(z)), np.ones(len(z))

    def test_the_trough_is_found(self, tmp_path):
        dcd, idx, m = self._bilayer(tmp_path, center=31.0)
        mid, com, off = bilayer_midplanes(dcd, idx, m)
        assert np.mod(mid, 60.0) == pytest.approx(31.0, abs=0.5)

    def test_a_bilayer_across_the_periodic_boundary(self, tmp_path):
        # center at z=2 in a 0..60 cell: half the lower leaflet is wrapped to the top
        dcd, idx, m = self._bilayer(tmp_path, center=2.0)
        mid, _, _ = bilayer_midplanes(dcd, idx, m)
        assert np.mod(mid[0] + 30.0, 60.0) - 30.0 == pytest.approx(2.0, abs=0.5)


def test_read_xst(tmp_path):
    p = tmp_path / 'r.xst'
    p.write_text('# NAMD extended system trajectory file\n'
                 '#$LABELS step a_x a_y a_z b_x b_y b_z c_x c_y c_z o_x o_y o_z s_x s_y s_z s_u s_v s_w\n'
                 '500 78 0 0 0 67 0 0 0 77.1 41.6 36.0 44.9 0 0 0 0 0 0\n')
    assert read_xst(str(p)) == {500: (77.1, 44.9)}
