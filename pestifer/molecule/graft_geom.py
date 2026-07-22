# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Geometry helpers for the graft presegment port: distance-based bond guessing and the pendant
(downstream) atom set of a rotatable bond -- the numpy analogue of VMD's autobond +
``PestiferCRot::rotate_pendant`` on a freshly loaded donor structure.

A graft donor is loaded from a raw PDB (no PSF connectivity), so -- like VMD's ``mol new`` -- bonds
are guessed from interatomic distances.  Glycan covalent bonds (~1.0-1.5 A) are well separated from
non-bonded contacts (>2.4 A), so the guess is unambiguous for the sugar chains a graft rotates.
"""
import numpy as np
from scipy.spatial import cKDTree

#: covalent radii (angstrom) for the elements that appear in glycan/protein donors
_COVRAD = {
    'H': 0.31, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'S': 1.05, 'P': 1.07,
    'F': 0.57, 'CL': 1.02, 'BR': 1.20, 'I': 1.39,
    'NA': 1.66, 'K': 2.03, 'CA': 1.76, 'MG': 1.41, 'ZN': 1.22,
}
#: added to the covalent-radii sum as the bond-length tolerance
_BOND_TOL = 0.45


def guess_bonds(coords, elems):
    """
    Guess covalent bonds from coordinates and element symbols.

    Two atoms are bonded when their distance is at most ``rcov_i + rcov_j + 0.45`` angstrom.

    Parameters
    ----------
    coords : numpy.ndarray
        ``(N, 3)`` atom coordinates.
    elems : sequence of str
        Per-atom element symbols.

    Returns
    -------
    list of set
        ``adj[i]`` is the set of atom indices bonded to atom ``i``.
    """
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    rad = np.array([_COVRAD.get(str(e).strip().upper(), 0.77) for e in elems])
    adj = [set() for _ in range(n)]
    if n < 2:
        return adj
    tree = cKDTree(coords)
    max_reach = 2.0 * float(rad.max()) + _BOND_TOL
    for i, j in tree.query_pairs(max_reach):
        if np.linalg.norm(coords[i] - coords[j]) <= rad[i] + rad[j] + _BOND_TOL:
            adj[i].add(j)
            adj[j].add(i)
    return adj


def pendant_indices(adj, i, j):
    """
    Atom indices of the pendant group downstream of the ``i -> j`` bond.

    Reproduces ``rotate_pendant``: start from ``j``'s neighbors (excluding ``i``) and grow by
    bond connectivity, never crossing ``i`` or ``j``.  ``j`` itself is excluded (it lies on the
    rotation axis, so it does not move regardless).

    Parameters
    ----------
    adj : list of set
        Bond adjacency (as from :func:`guess_bonds`).
    i, j : int
        The two bonded atom indices; the pendant is the ``j`` side of the cut.

    Returns
    -------
    set of int
        The pendant atom indices (empty if ``i`` and ``j`` are not bonded).
    """
    if j not in adj[i]:
        return set()
    pend = {k for k in adj[j] if k != i}
    frontier = list(pend)
    while frontier:
        nxt = []
        for a in frontier:
            for b in adj[a]:
                if b != i and b != j and b not in pend:
                    pend.add(b)
                    nxt.append(b)
        frontier = nxt
    return pend
