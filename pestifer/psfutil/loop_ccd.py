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


def dihedral_deg(p1, p2, p3, p4):
    """
    Signed dihedral (degrees, (-180, 180]) of the four points p1-p2-p3-p4 (each (3,)).

    Sign follows the IUPAC / VMD ``measure dihed`` convention (verified against VMD), so
    values are directly comparable to standard Ramachandran phi/psi (e.g. alpha ~ -63/-43).
    """
    b1 = np.asarray(p2, float) - np.asarray(p1, float)
    b2 = np.asarray(p3, float) - np.asarray(p2, float)
    b3 = np.asarray(p4, float) - np.asarray(p3, float)
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return float(np.degrees(np.arctan2(-y, x)))   # -y matches VMD's sign


def set_dihedral(coords, i, j, k, l, target_deg, moving_mask):
    """
    Rotate the atoms flagged by ``moving_mask`` about the j-k bond so that the dihedral
    ``coords[i]-coords[j]-coords[k]-coords[l]`` equals ``target_deg``.

    ``coords`` is (M, 3); ``i, j, k, l`` are indices into it; ``moving_mask`` is a length-M
    boolean array of the atoms downstream of the j-k bond (must include ``l``, must exclude
    ``i, j, k``). Returns a rotated copy; the input is not modified.
    """
    X = np.array(coords, dtype=float)
    current = dihedral_deg(X[i], X[j], X[k], X[l])
    delta = target_deg - current
    mask = np.asarray(moving_mask, dtype=bool)
    # With the VMD-convention dihedral above, a right-handed rotation of the downstream atoms
    # about the j->k axis increases the dihedral, so rotate by +delta about (X[k]-X[j]).
    X[mask] = rotate_points_about_axis(X[mask], X[j], X[k] - X[j], delta)
    return X

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


_BACKBONE = ('N', 'CA', 'C', 'O')


def backbone_from_pdb(pdb_path, chainID=None, segname=None):
    """
    Parse backbone atoms (N, CA, C, O) from a PDB into ``{resid: {atomname: xyz}}``.

    Selects records by ``segname`` (cols 73-76) when given, else by ``chainID`` (col 22),
    else all. resid is the integer residue sequence number (col 23-26). Later records win on
    a duplicate (resid, atomname) -- callers should pass an altloc-free / single-model PDB.
    """
    out = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            name = line[12:16].strip()
            if name not in _BACKBONE:
                continue
            if segname is not None and line[72:76].strip() != segname:
                continue
            if chainID is not None and line[21:22].strip() != chainID:
                continue
            try:
                resid = int(line[22:26])
                xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                continue
            out.setdefault(resid, {})[name] = xyz
    return out


def build_loop_ccd_problem(backbone, loop_resids, n_anchor_resid):
    """
    Assemble the CCD arrays for closing one loop, from a per-residue backbone map.

    Parameters
    ----------
    backbone : dict[int, dict[str, ndarray]]
        ``{resid: {atomname: xyz}}`` for at least the atoms named in ``N, CA, C, O`` of the
        loop residues plus the N-terminal anchor's C (needed as the phi reference of the first
        loop residue).
    loop_resids : sequence[int]
        The loop residues, N->C order (the movable residues).
    n_anchor_resid : int
        The resolved residue just before the loop (fixed); its C anchors the first phi.

    Returns
    -------
    dict with keys:
        coords       : (M, 3) backbone atoms, order = [n_anchor:C] + per loop residue N,CA,C,O
        row          : {(resid, atomname): row index into coords}
        bonds        : list of (a_row, b_row) rotatable backbone bonds, N->C
                       (each loop residue contributes phi = N-CA and psi = CA-C)
        moving_masks : list of bool arrays (len M), atoms downstream of each bond
        end_idx      : rows of the last loop residue's N, CA, C (the end effector)

    The caller supplies the closure ``target`` (the desired positions of ``end_idx``): the
    native positions for a delete-and-rebuild test, or the downstream anchor's geometry for a
    real build.
    """
    order = []          # (resid, atomname) in coords order
    coords = []
    coords.append(backbone[n_anchor_resid]['C']); order.append((n_anchor_resid, 'C'))
    for r in loop_resids:
        for an in _BACKBONE:
            if an in backbone[r]:
                coords.append(backbone[r][an]); order.append((r, an))
    coords = np.asarray(coords, dtype=float)
    row = {key: i for i, key in enumerate(order)}
    M = len(coords)

    def downstream_mask(after_row):
        m = np.zeros(M, dtype=bool)
        m[after_row + 1:] = True
        return m

    bonds, masks = [], []
    for r in loop_resids:
        n, ca, c = row[(r, 'N')], row[(r, 'CA')], row[(r, 'C')]
        # phi: rotate about N-CA -> moves everything after CA (C, O, downstream residues)
        bonds.append((n, ca)); masks.append(downstream_mask(ca))
        # psi: rotate about CA-C -> moves everything after C (O, downstream residues)
        bonds.append((ca, c)); masks.append(downstream_mask(c))
    last = loop_resids[-1]
    end_idx = [row[(last, 'N')], row[(last, 'CA')], row[(last, 'C')]]
    return {'coords': coords, 'row': row, 'bonds': bonds,
            'moving_masks': masks, 'end_idx': end_idx}


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
