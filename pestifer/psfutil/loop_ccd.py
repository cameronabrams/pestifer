# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Physics-based loop closure by cyclic coordinate descent (CCD).

Offline, deterministic, and free of VMD/NAMD: everything here operates on numpy
coordinate arrays. This is the geometry core of P1 of the loop-modeling redesign
(``docs/design/loop-modeling.md``): it closes a built (``guesscoord``) loop onto its
downstream anchor by distributing the closure across the loop's *own* backbone
dihedrals, rather than dragging one end across space (the steered-MD artifact this
replaces).

The core is deliberately protein-agnostic -- it moves index-addressed points about
index-addressed bond axes toward a fixed target -- so it can be unit-tested on synthetic
geometry. The task layer is responsible for mapping a protein loop's backbone atoms onto
these arrays.

Reference: A. A. Canutescu & R. L. Dunbrack Jr., "Cyclic coordinate descent: A robotics
algorithm for protein loop closure", Protein Science 12:963-972 (2003).
"""
import numpy as np

from ..util.coord import rotate_points_about_axis

# Analytic Ramachandran basins for the P1 initial guess: name -> (phi_deg, psi_deg, weight).
# Zero-data (no derived library yet -- that is P2). Centers are the conventional means of the
# major coil-accessible basins; weights are a coil-like mixture (helix/sheet/PPII dominant,
# left-handed rare).
RAMACHANDRAN_BASINS = {
    'alphaR': (-63.0, -43.0, 0.42),   # right-handed alpha
    'beta':   (-135.0, 135.0, 0.33),  # beta / extended
    'ppII':   (-75.0, 145.0, 0.20),   # polyproline II
    'alphaL': (60.0, 45.0, 0.05),     # left-handed alpha (mostly Gly)
}


def sample_backbone_dihedrals(rng, n, jitter_deg=10.0):
    """
    Draw ``n`` (phi, psi) pairs, in degrees, from the :data:`RAMACHANDRAN_BASINS` mixture.

    Deterministic for a fixed ``rng`` (a ``numpy.random.Generator``): a build seeded the
    same way samples the same starting backbone. Each residue picks a basin by weight, then
    a Gaussian jitter (std ``jitter_deg``) is added to each angle. Returns an ``(n, 2)``
    array of ``[phi, psi]`` in degrees (not wrapped; callers that need [-180, 180) should
    wrap).
    """
    names = list(RAMACHANDRAN_BASINS)
    weights = np.array([RAMACHANDRAN_BASINS[k][2] for k in names], dtype=float)
    weights /= weights.sum()
    picks = rng.choice(len(names), size=n, p=weights)
    centers = np.array([[RAMACHANDRAN_BASINS[names[p]][0],
                         RAMACHANDRAN_BASINS[names[p]][1]] for p in picks], dtype=float)
    return centers + rng.normal(0.0, jitter_deg, size=(n, 2))


def optimal_ccd_angle(moving, target, pivot, axis, weights=None):
    """
    Closed-form rotation angle (degrees, right-handed about ``axis``) that minimizes
    ``sum_k w_k |R(theta) moving_k - target_k|^2`` for a rotation about the line through
    ``pivot`` along ``axis``.

    Because each moving atom's component perpendicular to the axis rotates in-plane, the
    objective's theta-dependence is ``A cos(theta) + B sin(theta)`` (up to constants), whose
    maximizer is ``theta* = atan2(B, A)`` with ``A = sum_k w_k (r_k . f_k)`` and
    ``B = sum_k w_k ((k_hat x r_k) . f_k)``, where ``r_k`` / ``f_k`` are the perpendicular
    components (relative to the axis) of ``moving_k`` / ``target_k`` about ``pivot``.

    Parameters
    ----------
    moving, target : array_like, shape (K, 3)
        Current and desired positions of the K end-effector atoms.
    pivot : array_like, shape (3,)
    axis : array_like, shape (3,)  (need not be unit length; a zero axis yields 0.0)
    weights : array_like, shape (K,), optional

    Returns
    -------
    float
        The optimal angle in degrees.
    """
    moving = np.asarray(moving, dtype=float)
    target = np.asarray(target, dtype=float)
    pivot = np.asarray(pivot, dtype=float)
    axis = np.asarray(axis, dtype=float)
    n = np.linalg.norm(axis)
    if n == 0.0 or len(moving) == 0:
        return 0.0
    k = axis / n
    w = np.ones(len(moving)) if weights is None else np.asarray(weights, dtype=float)
    r = moving - pivot
    r = r - np.outer(r @ k, k)          # perpendicular component of moving atoms
    f = target - pivot
    f = f - np.outer(f @ k, k)          # perpendicular component of targets
    A = float(np.sum(w * np.sum(r * f, axis=1)))
    kxr = np.cross(np.broadcast_to(k, r.shape), r)
    B = float(np.sum(w * np.sum(kxr * f, axis=1)))
    return float(np.degrees(np.arctan2(B, A)))


def end_rmsd(coords, end_idx, target):
    """RMSD (Angstrom) between the end-effector atoms ``coords[end_idx]`` and ``target``."""
    d = np.asarray(coords, dtype=float)[end_idx] - np.asarray(target, dtype=float)
    return float(np.sqrt(np.mean(np.sum(d * d, axis=1))))


def ccd_close(coords, bonds, moving_masks, end_idx, target,
              weights=None, max_iters=500, tol=0.08, stall_tol=1e-6):
    """
    Close a loop onto a fixed target by cyclic coordinate descent.

    Sweeps the rotatable ``bonds`` in order (N->C is the natural order for a loop); for each
    bond it applies :func:`optimal_ccd_angle` to bring the end-effector atoms onto ``target``,
    rotating only the atoms flagged by that bond's mask. Repeats until the end RMSD falls
    below ``tol`` or progress stalls.

    Parameters
    ----------
    coords : array_like, shape (M, 3)
        Loop atom coordinates (not modified; a copy is returned).
    bonds : list[tuple[int, int]]
        ``(a, b)`` index pairs into ``coords``; each defines a rotation axis (the line
        through ``coords[a]`` along ``coords[b] - coords[a]``).
    moving_masks : list[array_like[bool]], each length M
        For each bond, which atoms rotate with it. Every end-effector atom must be movable by
        every bond it is downstream of; typically a bond's mask is "all atoms C-terminal to
        it", which includes ``end_idx``.
    end_idx : array_like[int]
        Indices of the end-effector atoms whose landing on ``target`` defines closure.
    target : array_like, shape (len(end_idx), 3)
    weights : array_like, optional
        Per-end-atom weights for the least-squares objective.
    max_iters : int
    tol : float
        Convergence RMSD in Angstrom.
    stall_tol : float
        Stop early if an iteration improves the RMSD by less than this.

    Returns
    -------
    (coords_closed, final_rmsd, iterations) : (np.ndarray, float, int)
    """
    X = np.array(coords, dtype=float)
    target = np.asarray(target, dtype=float)
    end_idx = np.asarray(end_idx, dtype=int)
    last = end_rmsd(X, end_idx, target)
    for it in range(1, max_iters + 1):
        for (a, b), mask in zip(bonds, moving_masks):
            pivot = X[a]
            axis = X[b] - X[a]
            theta = optimal_ccd_angle(X[end_idx], target, pivot, axis, weights)
            mask = np.asarray(mask, dtype=bool)
            X[mask] = rotate_points_about_axis(X[mask], pivot, axis, theta)
        r = end_rmsd(X, end_idx, target)
        if r < tol or (last - r) < stall_tol:
            return X, r, it
        last = r
    return X, end_rmsd(X, end_idx, target), max_iters
