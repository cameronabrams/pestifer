# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""Knots in a protein backbone: the one loop-modeling failure no minimizer or MD run can undo.

A modeled loop that passes one strand through another (or through itself) is reported at closure
time by a proxy -- two non-adjacent alpha carbons closer than 3 A (``THREADED`` in
``LoopClosureMixin``) -- which is a close approach, not proof that anything crossed.  What a real
strand passage leaves behind, once minimization has pushed the strands apart on whichever side
they happened to be, is a change of crossing; and a crossing change in a chain is exactly what
ties a knot.  Natively knotted proteins exist but are rare (well under 1% of the PDB) and known,
so a knot in a built chain is close to proof of a modeling error.

The test, on each covalently continuous alpha-carbon trace:

1. Close the open chain far outside the protein: extend each terminus radially from the chain's
   centroid to a distance many times the chain's size, and join the two far points through a
   third far point.  This is the usual "closure to infinity"; it is deterministic, and for a
   terminus buried in the protein it can itself tie or untie a knot, which is the known
   ambiguity of knotting in open chains (the same one KnotProt resolves by sampling closures).
2. Simplify with Koniaris-Muthukumar-Taylor reduction: delete a vertex whenever the triangle it
   makes with its two neighbors is pierced by no other edge.  Each deletion is an isotopy, so the
   knot type is unchanged; an unknotted chain usually collapses to a triangle.
3. Compute the knot determinant ``|Delta(-1)|`` from a generic projection of what is left: the
   Fox-coloring matrix (row per crossing: 2 at the over-arc, -1 at each under-arc, independent of
   crossing sign), any one minor, in exact integer arithmetic.  The unknot gives 1, the trefoil 3,
   the figure-eight 5.  A few knots also give 1 and are missed; none is a plausible modeling
   artifact.

Only covalent backbone continuity is used: disulfides are not followed, so a loop threaded
through a cystine-closed lasso is not seen here, and two separate chains entangled with each
other have no knot invariant to compute.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

_EPS = 1e-9


@dataclass(frozen=True)
class Knot:
    segname: str
    first: str            # resid of the chain's first and last alpha carbon
    last: str
    determinant: int
    kind: str = 'knot'

    @property
    def resid(self):
        return self.first

    @property
    def key(self):
        return ('knot', self.segname, self.first)

    def describe(self):
        name = {3: 'a trefoil', 5: 'a figure-eight'}.get(self.determinant, 'a knot')
        return (f'{self.segname}:{self.first}-{self.last} backbone is tied in {name} '
                f'(knot determinant {self.determinant}) -- a strand passed through another; '
                f'no minimization can undo it')


def backbone_chains(peptides):
    """Covalently continuous alpha-carbon runs, as lists of atom rows in chain order.

    ``peptides`` are :class:`~pestifer.psfutil.stereo.StereoChecker` peptide tuples, whose
    fourth and last fields are the alpha-carbon rows of the two residues a peptide bond joins.
    """
    nxt, prv = {}, {}
    for t in peptides:
        nxt[t[3]] = t[6]
        prv[t[6]] = t[3]
    chains = []
    for start in nxt:
        if start in prv:
            continue
        run = [start]
        while run[-1] in nxt and len(run) <= len(nxt):
            run.append(nxt[run[-1]])
        chains.append(run)
    return chains


def close_to_infinity(P):
    """The open polyline ``P`` closed through three points far outside it (see module doc)."""
    c = P.mean(axis=0)
    R = 100.0 * (np.linalg.norm(P - c, axis=1).max() + 1.0)

    def out(p):
        d = p - c
        n = np.linalg.norm(d)
        return d / n if n > _EPS else np.array([1.0, 0.0, 0.0])
    d0, d1 = out(P[0]), out(P[-1])
    mid = d0 + d1
    if np.linalg.norm(mid) < 1e-3:          # termini leave in opposite directions
        mid = np.cross(d0, [1.0, 0.0, 0.0] if abs(d0[0]) < 0.9 else [0.0, 1.0, 0.0])
    mid /= np.linalg.norm(mid)
    return np.vstack([P, c + R * d1, c + 2 * R * mid, c + R * d0])


def _triangle_hits(a, b, c, S0, S1):
    """Which segments ``S0[k]->S1[k]`` pass through triangle ``abc`` (Moller-Trumbore)."""
    e1, e2 = b - a, c - a
    D = S1 - S0
    h = np.cross(D, e2)
    det = h @ e1
    ok = np.abs(det) > _EPS
    inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
    s = S0 - a
    u = inv * np.einsum('ij,ij->i', s, h)
    q = np.cross(s, e1)
    v = inv * np.einsum('ij,ij->i', D, q)
    t = inv * (q @ e2)
    return ok & (u >= -_EPS) & (v >= -_EPS) & (u + v <= 1 + _EPS) & (t >= -_EPS) & (t <= 1 + _EPS)


def kmt_reduce(Q):
    """Koniaris-Muthukumar-Taylor reduction of the closed polygon ``Q`` (vertices in order)."""
    Q = np.array(Q, dtype=float)
    changed = True
    while changed and len(Q) > 3:          # passes until one deletes nothing
        changed = False
        i = 0
        while i < len(Q) and len(Q) > 3:
            n = len(Q)
            # edge k runs Q[k] -> Q[k+1]; skip the four that touch the triangle's vertices
            keep = np.ones(n, dtype=bool)
            keep[[(i - 2) % n, (i - 1) % n, i, (i + 1) % n]] = False
            if keep.any() and _triangle_hits(Q[(i - 1) % n], Q[i], Q[(i + 1) % n],
                                             Q[keep], np.roll(Q, -1, axis=0)[keep]).any():
                i += 1
                continue
            Q = np.delete(Q, i, axis=0)
            changed = True
    return Q


def _det_int(M):
    """Exact determinant of an integer matrix (Bareiss fraction-free elimination)."""
    M = [list(map(int, row)) for row in M]
    n = len(M)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for k in range(n - 1):
        if M[k][k] == 0:
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    M[k], M[r] = M[r], M[k]
                    sign = -sign
                    break
            else:
                return 0
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def knot_determinant(Q, seed=1):
    """``|Delta(-1)|`` of the closed polygon ``Q`` (1 for the unknot)."""
    Q = np.asarray(Q, dtype=float)
    n = len(Q)
    if n < 6:                                   # a knot needs at least six edges
        return 1
    rng = np.random.default_rng(seed)
    z = rng.normal(size=3)
    z /= np.linalg.norm(z)
    x = np.cross(z, [1.0, 0.0, 0.0] if abs(z[0]) < 0.9 else [0.0, 1.0, 0.0])
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    P2 = np.stack([Q @ x, Q @ y], axis=1)
    H = Q @ z
    crossings = []                              # (pos_over, pos_under)
    for i in range(n):
        p, r = P2[i], P2[(i + 1) % n] - P2[i]
        for j in range(i + 2, n):
            if (j + 1) % n == i:
                continue
            q, s = P2[j], P2[(j + 1) % n] - P2[j]
            den = r[0] * s[1] - r[1] * s[0]
            if abs(den) < _EPS:
                continue
            w = q - p
            t = (w[0] * s[1] - w[1] * s[0]) / den
            u = (w[0] * r[1] - w[1] * r[0]) / den
            if not (0.0 < t < 1.0 and 0.0 < u < 1.0):
                continue
            hi = H[i] + t * (H[(i + 1) % n] - H[i])
            hj = H[j] + u * (H[(j + 1) % n] - H[j])
            crossings.append((i + t, j + u) if hi > hj else (j + u, i + t))
    N = len(crossings)
    if N == 0:
        return 1
    unders = sorted(c[1] for c in crossings)

    def arc(pos):                               # arc k runs from under-crossing k-1 to k
        return int(np.searchsorted(unders, pos)) % N
    M = np.zeros((N, N), dtype=np.int64)
    for k, (over, under) in enumerate(crossings):
        u = int(np.searchsorted(unders, under))
        M[k, arc(over)] += 2
        M[k, u % N] -= 1
        M[k, (u + 1) % N] -= 1
    return abs(_det_int(M[:-1, :-1].tolist()))


def find_knots(xyz, chains, label) -> list[Knot]:
    """Every knotted chain among ``chains`` (lists of alpha-carbon rows) in frame ``xyz``.

    ``label(row)`` gives ``(segname, resid)`` for a row.
    """
    out = []
    for run in chains:
        if len(run) < 6:
            continue
        Q = kmt_reduce(close_to_infinity(xyz[run]))
        det = knot_determinant(Q)
        if det != 1:
            seg, first = label(run[0])
            out.append(Knot(seg, first, label(run[-1])[1], det))
    return out
