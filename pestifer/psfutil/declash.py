# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Numpy declashers -- the VMD-free replacement for ``PestiferDeclash::declash_pendant`` and
``declash_loop`` (declash.tcl).

Declashing is *pre-conditioning*: it rotates rotatable bonds (pendant sub-branches) or backbone
phi angles (gap loops) to pull model-built / grafted atoms out of steric overlap with the rest of
the structure before minimization, which then relaxes the result.  It is not structure
prediction, so the exact poses do not matter -- only that the destabilizing heavy-atom clash
count drops.  Both routines are the same simple greedy descent as the Tcl (accept a random move
only if it strictly reduces the clash count), scored with a ``scipy.spatial.cKDTree`` heavy-atom
contact counter, and driven by a caller-supplied seeded RNG so declashed builds are reproducible
(VMD's ``rand()`` was unseeded).
"""
import numpy as np
from scipy.spatial import cKDTree

from ..util.coord import rotate_points_about_axis

#: discrete rotation increments (degrees) tried for a pendant bond, matching declash.tcl
_PENDANT_DEGS = np.array([-120.0, -60.0, 60.0, 120.0, 180.0])


def _contact_count(env_tree: cKDTree, pts: np.ndarray, clashdist: float) -> int:
    """Number of (query point, env atom) pairs within ``clashdist`` -- the clash score."""
    if len(pts) == 0:
        return 0
    return int(env_tree.query_ball_point(pts, clashdist, return_length=True).sum())


def declash_pendant(coords, pendant_heavy, env_heavy, bonds, movers, maxcycles, clashdist, rng):
    """
    Greedily rotate a pendant group's rotatable bonds to reduce heavy-atom clashes with the rest
    of the structure.  ``coords`` is mutated in place.

    Parameters
    ----------
    coords : numpy.ndarray
        ``(N, 3)`` coordinates of the whole system.
    pendant_heavy : numpy.ndarray
        Row indices of the pendant's heavy atoms (the atoms whose clashes we score).
    env_heavy : numpy.ndarray
        Row indices of every *other* heavy atom (the environment; static during this call).
    bonds : list of (int, int)
        Rotatable bonds as pairs of atom row indices; the axis runs from the first to the second.
    movers : list of numpy.ndarray
        Per-bond arrays of atom row indices that rotate about that bond.
    maxcycles : int
        Maximum Monte-Carlo cycles (one bond rotation each).
    clashdist : float
        Heavy-atom clash cutoff (angstrom).
    rng : numpy.random.Generator
        Seeded RNG driving the random bond/angle choices.

    Returns
    -------
    int
        The final clash count (0 if fully resolved).
    """
    if not bonds or len(pendant_heavy) == 0 or len(env_heavy) == 0:
        return 0
    env_tree = cKDTree(coords[env_heavy])
    ncontacts = _contact_count(env_tree, coords[pendant_heavy], clashdist)
    if ncontacts == 0:
        return 0
    nb = len(bonds)
    for _ in range(maxcycles):
        ridx = int(rng.integers(nb))
        bi, bj = bonds[ridx]
        mv = movers[ridx]
        deg = float(rng.choice(_PENDANT_DEGS))
        pivot, axis = coords[bi], coords[bj] - coords[bi]
        coords[mv] = rotate_points_about_axis(coords[mv], pivot, axis, deg)
        new = _contact_count(env_tree, coords[pendant_heavy], clashdist)
        if new >= ncontacts:                       # reject: rotate back
            coords[mv] = rotate_points_about_axis(coords[mv], pivot, axis, -deg)
        else:
            ncontacts = new
            if ncontacts == 0:
                break
    return ncontacts


def declash_loop(coords, residues, env_heavy, maxcycles, clashdist, rng, jitter=120.0):
    """
    Greedily wiggle a model-built gap loop's backbone phi angles to reduce clashes with the
    surrounding structure.  ``coords`` is mutated in place.

    Parameters
    ----------
    coords : numpy.ndarray
        ``(N, 3)`` coordinates of the whole system.
    residues : list of dict
        One entry per loop residue, innermost first, each with keys ``pivot`` (N atom row),
        ``axis_to`` (CA atom row), ``movers`` (row indices rotated by this residue's phi -- this
        residue's non-N/HN atoms plus every downstream loop atom), and ``frag_heavy`` (heavy row
        indices of this residue through the loop C-terminus -- the atoms whose clashes are scored).
    env_heavy : numpy.ndarray
        Row indices of the non-loop heavy atoms (the environment; static during this call).
    maxcycles : int
        Maximum Monte-Carlo cycles per residue.
    clashdist : float
        Heavy-atom clash cutoff (angstrom).
    rng : numpy.random.Generator
        Seeded RNG.
    jitter : float
        Half-width (degrees) of the uniform phi perturbation.

    Returns
    -------
    int
        The final total clash count summed over the loop residues processed.
    """
    if len(env_heavy) == 0:
        return 0
    env_tree = cKDTree(coords[env_heavy])
    total = 0
    for res in residues:
        frag = res['frag_heavy']
        movers = res['movers']
        con = _contact_count(env_tree, coords[frag], clashdist)
        for _ in range(maxcycles):
            if con == 0:
                break
            deg = float((1.0 - 2.0 * rng.random()) * jitter)
            pivot, axis = coords[res['pivot']], coords[res['axis_to']] - coords[res['pivot']]
            save = coords[movers].copy()
            coords[movers] = rotate_points_about_axis(coords[movers], pivot, axis, deg)
            new = _contact_count(env_tree, coords[frag], clashdist)
            if new < con:
                con = new
            else:
                coords[movers] = save                # reject
        total += con
    return total
