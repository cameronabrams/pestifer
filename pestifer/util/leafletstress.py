# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Per-leaflet surface tension from a pressure-profile trajectory, split at the bilayer midplane.

Two things here replace assumptions the differential-stress diagnostic used to make:

**Where the midplane is.**  It used to be "the central slab".  Nothing keeps a bilayer centered in
its cell, and the leaflet tensions are a split of a profile whose integrand is large in the
headgroup and tail regions, so one slab of misplacement (~2 A at 40 slabs) moved ex17's Dgamma by
8-15 mN/m -- more than the signal.  Here the midplane is measured from the coordinates, per frame:
the lipid center of mass tracks drift exactly, and a fixed offset from it to the methyl trough
(the lipid-density minimum between the leaflets, located once from all frames pooled) puts the
split where the leaflets meet, not where their mass balances -- those differ in an asymmetric
bilayer.

**Where each leaflet ends.**  The lower leaflet is integrated from the midplane down half a cell
and the upper from the midplane up half a cell, periodically, so the outer boundaries sit in the
middle of the water whichever way the bilayer has drifted.

NAMD's slab geometry, from the 3.0.3 source (``ComputeNonbondedPair.C:219``, ``ComputeEwald.C:135``,
the same for both profile kinds): slab 0 starts at ``origin.z - c/2`` and slabs are ``c/nslabs``
thick.  The origin is not in a DCD; it comes from the ``.xst``.
"""
from __future__ import annotations

import numpy as np

from .dcd import iter_frames
from .densityprofile import _parse_psf, classify_species

#: bin width (A) of the pooled lipid-density histogram the methyl trough is located in
TROUGH_BIN = 0.25
#: the trough is searched for within this distance (A) of the lipid center of mass
TROUGH_WINDOW = 8.0


def read_xst(path: str) -> dict[int, tuple[float, float]]:
    """``{step: (c_z, o_z)}`` from a NAMD ``.xst``: the cell depth and the origin's z, in A."""
    out = {}
    with open(path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            v = line.split()
            out[int(float(v[0]))] = (float(v[9]), float(v[12]))
    return out


def lipid_atoms(psf: str) -> tuple[np.ndarray, np.ndarray]:
    """Indices and masses of the lipid HEAVY atoms of ``psf``, in PSF order."""
    resnames, masses = _parse_psf(psf)
    sel = (classify_species(resnames) == 'lipid') & (masses > 1.5)
    idx = np.flatnonzero(sel)
    if len(idx) == 0:
        raise ValueError(f'{psf}: no lipid atoms; cannot locate a bilayer midplane')
    return idx, masses[idx]


def _unwrapped_about_center(z: np.ndarray, c: float) -> tuple[float, np.ndarray]:
    """Circular-mean center of ``z`` in a period-``c`` box, and ``z`` unwrapped about it."""
    th = 2.0 * np.pi * z / c
    center = np.arctan2(np.sin(th).mean(), np.cos(th).mean()) * c / (2.0 * np.pi)
    rel = z - center
    rel -= c * np.round(rel / c)
    return center, rel


def bilayer_midplanes(dcd: str, idx: np.ndarray, masses: np.ndarray, frames=None):
    """Per-frame midplane z (A, in the DCD's frame) and the trough offset from the lipid COM.

    Returns ``(midplane, com, offset)``: ``midplane = com + offset``.
    """
    coms, rels = [], []
    for _, cell, xyz in iter_frames(dcd, atoms=idx, frames=frames):
        if cell is None:
            raise ValueError(f'{dcd}: no unit cell in the DCD (write it with DCDUnitCell yes)')
        center, rel = _unwrapped_about_center(xyz[:, 2], cell[2])
        shift = float(np.average(rel, weights=masses))
        coms.append(center + shift)
        rels.append(rel - shift)
    coms = np.asarray(coms)
    pooled = np.concatenate(rels)
    w = np.tile(masses, len(rels))
    edges = np.arange(-TROUGH_WINDOW, TROUGH_WINDOW + TROUGH_BIN / 2, TROUGH_BIN)
    hist, _ = np.histogram(pooled, bins=edges, weights=w)
    k = 4                                    # ~1 A boxcar: the trough, not single-bin noise
    smooth = np.convolve(hist, np.ones(k) / k, mode='same')
    centers = 0.5 * (edges[:-1] + edges[1:])
    inner = slice(k, len(smooth) - k)        # avoid the boxcar's edge artifacts
    j = int(np.argmin(smooth[inner])) + k
    # The trough is the contiguous run of bins around the minimum that sit within 10% (of the
    # window's peak-to-minimum range) of it, located at its depth-weighted centroid.  A bare argmin
    # takes the FIRST bin of a flat-bottomed trough, which put a synthetic 4 A gap's midplane at its
    # edge, 1.2 A off.
    lo_v, hi_v = smooth[j], smooth[inner].max()
    thr = lo_v + 0.1 * (hi_v - lo_v)
    a = j
    while a - 1 >= k and smooth[a - 1] <= thr:
        a -= 1
    b = j
    while b + 1 < len(smooth) - k and smooth[b + 1] <= thr:
        b += 1
    w = thr - smooth[a:b + 1] + 1e-12
    offset = float(np.average(centers[a:b + 1], weights=w))
    return coms + offset, coms, offset


def slab_tensions(profiles: np.ndarray, c_z: np.ndarray) -> np.ndarray:
    """``(nframes, nslabs)`` slab tensions, mN/m: ``0.01 * (Pzz - (Pxx+Pyy)/2) * dz``.

    ``profiles`` is ``(nframes, nslabs, 3)`` in bar; ``c_z`` the per-frame cell depth in A.
    """
    p = np.asarray(profiles, dtype=float)
    dz = np.asarray(c_z, dtype=float) / p.shape[1]
    return 0.01 * (p[..., 2] - 0.5 * (p[..., 0] + p[..., 1])) * dz[:, None]


def _periodic_cumulative(g: np.ndarray, x: float) -> float:
    """Integral of the periodic step function ``g`` (one value per unit slab) from 0 to ``x``."""
    s = len(g)
    whole, frac = divmod(x, 1.0)
    whole = int(whole)
    n, r = divmod(whole, s)
    total = n * g.sum() + g[:r].sum()
    return float(total + frac * g[r % s])


def leaflet_tensions(gamma_slab: np.ndarray, midplane: np.ndarray, origin_z: np.ndarray,
                     c_z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-frame ``(lower, upper)`` leaflet tensions, each integrated over half a cell.

    ``gamma_slab`` is ``(nframes, nslabs)`` from :func:`slab_tensions`; the other arguments are
    per frame, in A.  The split point is the midplane converted to slab units, fractional.
    """
    g = np.asarray(gamma_slab, dtype=float)
    nf, s = g.shape
    lower, upper = np.empty(nf), np.empty(nf)
    for f in range(nf):
        dz = c_z[f] / s
        u = ((midplane[f] - (origin_z[f] - 0.5 * c_z[f])) / dz) % s + s   # keep u - s/2 > 0
        F = lambda x: _periodic_cumulative(g[f], x)
        lower[f] = F(u) - F(u - s / 2.0)
        upper[f] = F(u + s / 2.0) - F(u)
    return lower, upper
