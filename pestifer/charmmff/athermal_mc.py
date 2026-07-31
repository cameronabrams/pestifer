# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Athermal (repulsive-only) Metropolis Monte Carlo sampler for single-lipid tail conformers.

The grid membrane packer needs a *fluid-like* ensemble of single-molecule conformers: shapes
with the moderately splayed, gauche-rich tails a lipid holds in a fluid bilayer.  Vacuum
single-molecule MD cannot produce these -- attractive van der Waals forces drive the isolated
chains into a collapsed globule, and stretching them to avoid that overshoots into thin
all-trans rods.  Neither packs to fluid density; both are the wrong ensemble (see
``docs/design/lipid-conformer-generation.md``).

This module breaks that chicken-and-egg with **athermal MC**: a purely excluded-volume
(hard-sphere) potential with *no* attractive term.  Removing the attraction removes the
collapse driving force, so the tails explore torsional entropy freely and settle into a
melt-like distribution -- a good model for the packing-dominated hydrophobic core.  "Athermal"
means geometry-only (the T->infinity limit): a move is accepted iff it creates no hard-sphere
overlap and keeps the confined atoms inside a lateral **cylinder**, so we sample configurational
entropy subject to excluded volume, exactly the missing quantity.

The engine here is deliberately force-field agnostic: it works on coordinates, per-atom
hard-sphere radii, a precomputed list of rotatable bonds (each with the tip-side atom set it
pivots), a pair-exclusion set, and a confinement cylinder.  The CHARMM-specific wiring that
builds those inputs for a real lipid lives elsewhere; this keeps the MC core small and unit
testable on a synthetic chain.

Moves are **dihedral pivots**: a rotatable bond's tip-side subtree is rigidly rotated about the
bond axis.  Bonds and angles are therefore preserved exactly (only torsions change), so every
sampled conformer retains the input CHARMM geometry to machine precision.
"""
from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

logger = logging.getLogger(__name__)

# Rmin (the CHARMM pair-potential minimum separation) is 2**(1/6) * sigma, where sigma is the
# Lennard-Jones/WCA collision diameter.  A hard sphere should exclude at ~sigma, not at the
# (fatter) energy minimum, so per-atom hard radii default to Rmin/2 scaled down to sigma/2.
_SIGMA_FROM_RMIN = 2.0 ** (-1.0 / 6.0)


@dataclass
class RotatableBond:
    """A single dihedral-pivot degree of freedom.

    Rotating about the directed axis ``coords[b] - coords[a]`` rigidly moves ``moving`` (the
    tip-side subtree, i.e. every atom that becomes disconnected from ``a`` when the ``a-b`` bond
    is cut), leaving ``a`` and everything on its side fixed.  Because the rotation is rigid, all
    bond lengths and bond angles are preserved; only torsions about ``a-b`` change.
    """
    a: int
    """Index of the pivot's anchor-side atom (stays fixed)."""
    b: int
    """Index of the pivot's tip-side atom (on the rotating axis, moves with ``moving``)."""
    moving: np.ndarray
    """Integer indices of the atoms rigidly rotated (the tip-side subtree, including ``b``)."""


@dataclass
class MoleculeMC:
    """Everything the athermal MC engine needs about one molecule.

    Parameters
    ----------
    coords : np.ndarray
        ``(N, 3)`` starting coordinates (typically an extended, valid-geometry conformer).
    radii : np.ndarray
        ``(N,)`` per-atom hard-sphere radii.  A pair ``(i, j)`` overlaps when their separation
        is below ``radii[i] + radii[j]``.
    rotatable : list[RotatableBond]
        The dihedral-pivot degrees of freedom (typically the acyl-tail C-C bonds).
    exclusions : set[frozenset[int]]
        Atom-index pairs to ignore in overlap tests -- the 1-2 (bonded) and 1-3 (shared-angle)
        neighbors, which are always within contact by construction.
    axis_point, axis_dir : np.ndarray
        A point on, and the unit direction of, the confinement cylinder's axis (the molecule's
        head->tail reference axis).
    cylinder_radius : float
        Radius of the confinement cylinder.  Confined atoms whose radial distance from the axis
        exceeds this are rejected.  ``inf`` disables confinement.
    confined : np.ndarray
        ``(N,)`` boolean mask of atoms subject to the cylinder (typically the tail atoms).
    """
    coords: np.ndarray
    radii: np.ndarray
    rotatable: list[RotatableBond]
    exclusions: set = field(default_factory=set)
    axis_point: np.ndarray = None
    axis_dir: np.ndarray = None
    cylinder_radius: float = float('inf')
    confined: np.ndarray = None

    def __post_init__(self):
        self.coords = np.ascontiguousarray(self.coords, dtype=float)
        self.radii = np.asarray(self.radii, dtype=float)
        n = len(self.coords)
        if self.axis_point is None:
            self.axis_point = np.zeros(3)
        else:
            self.axis_point = np.asarray(self.axis_point, dtype=float)
        if self.axis_dir is None:
            self.axis_dir = np.array([0.0, 0.0, 1.0])
        else:
            d = np.asarray(self.axis_dir, dtype=float)
            self.axis_dir = d / np.linalg.norm(d)
        if self.confined is None:
            self.confined = np.ones(n, dtype=bool)
        else:
            self.confined = np.asarray(self.confined, dtype=bool)


def _rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """Rodrigues rotation matrix for angle ``theta`` about unit ``axis``."""
    x, y, z = axis
    c, s = np.cos(theta), np.sin(theta)
    C = 1.0 - c
    return np.array([
        [c + x * x * C,     x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, c + y * y * C,     y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
    ])


def radial_distances(coords: np.ndarray, point: np.ndarray, direction: np.ndarray) -> np.ndarray:
    """Perpendicular distance of each row of ``coords`` from the line ``point + t*direction``."""
    rel = coords - point
    axial = rel @ direction
    perp = rel - np.outer(axial, direction)
    return np.linalg.norm(perp, axis=1)


def _apply_pivot(coords: np.ndarray, bond: RotatableBond, theta: float) -> np.ndarray:
    """Return a copy of ``coords`` with ``bond.moving`` rotated by ``theta`` about the bond axis."""
    axis = coords[bond.b] - coords[bond.a]
    norm = np.linalg.norm(axis)
    if norm < 1e-9:
        return coords.copy()
    R = _rotation_matrix(axis / norm, theta)
    new = coords.copy()
    pivot = coords[bond.a]
    new[bond.moving] = (coords[bond.moving] - pivot) @ R.T + pivot
    return new


def _has_overlap(coords: np.ndarray, radii: np.ndarray, moving: np.ndarray,
                 exclusions: set) -> bool:
    """True if any non-excluded pair involving a moved atom is within hard-sphere contact.

    Only pairs with at least one moved atom can newly clash: the moving subtree is rigid, so
    moved-moved separations are unchanged from the (valid) pre-move state, and stationary-
    stationary pairs did not move at all.  We therefore query the full structure once and keep
    only candidate pairs that touch the moved set.
    """
    rmax = radii.max()
    tree = cKDTree(coords)
    moved = set(int(i) for i in moving)
    for i, j in tree.query_pairs(2.0 * rmax):
        if i not in moved and j not in moved:
            continue
        if frozenset((i, j)) in exclusions:
            continue
        if np.linalg.norm(coords[i] - coords[j]) < radii[i] + radii[j]:
            return True
    return False


def run_mc(mol: MoleculeMC, nsamples: int = 10, n_equil: int = 2000, n_decorr: int = 200,
           max_angle: float = np.pi, seed: int = None) -> list[np.ndarray]:
    """Run athermal dihedral-pivot MC and return ``nsamples`` decorrelated conformers.

    Starting from ``mol.coords``, propose pivots on random rotatable bonds by random angles in
    ``[-max_angle, max_angle]``; accept iff the result is overlap-free and keeps confined atoms
    inside the cylinder (the athermal, T->infinity criterion -- energy never enters).  After
    ``n_equil`` proposals to melt the (typically rod-like) starting conformer, collect one
    conformer every ``n_decorr`` proposals.

    Parameters
    ----------
    mol : MoleculeMC
        The molecule and its confinement.  Not mutated.
    nsamples : int
        Number of conformers to return.
    n_equil : int
        Proposals to run (and discard) before the first sample, to decorrelate from the start.
    n_decorr : int
        Proposals between retained samples.
    max_angle : float
        Half-width (radians) of the uniform pivot-angle proposal.  ``pi`` (default) proposes the
        full torsional circle for aggressive mixing.
    seed : int, optional
        Seed for the sampler's RNG (reproducible output).

    Returns
    -------
    list[np.ndarray]
        ``nsamples`` coordinate arrays, each ``(N, 3)``.
    """
    if not mol.rotatable:
        logger.warning('run_mc: molecule has no rotatable bonds; returning copies of the input')
        return [mol.coords.copy() for _ in range(nsamples)]

    rng = np.random.default_rng(seed)
    coords = mol.coords.copy()
    nbonds = len(mol.rotatable)
    samples: list[np.ndarray] = []
    n_accept = n_attempt = 0

    total_proposals = n_equil + nsamples * n_decorr
    for step in range(1, total_proposals + 1):
        bond = mol.rotatable[rng.integers(nbonds)]
        theta = rng.uniform(-max_angle, max_angle)
        trial = _apply_pivot(coords, bond, theta)
        n_attempt += 1
        # confinement: only the moved, confined atoms can newly leave the cylinder
        conf_moved = bond.moving[mol.confined[bond.moving]]
        if len(conf_moved) and np.isfinite(mol.cylinder_radius):
            if radial_distances(trial[conf_moved], mol.axis_point,
                                 mol.axis_dir).max() > mol.cylinder_radius:
                continue
        if _has_overlap(trial, mol.radii, bond.moving, mol.exclusions):
            continue
        coords = trial
        n_accept += 1
        if step > n_equil and (step - n_equil) % n_decorr == 0:
            samples.append(coords.copy())

    logger.debug(f'run_mc: {n_accept}/{n_attempt} moves accepted '
                 f'({100.0 * n_accept / max(n_attempt, 1):.1f}%); collected {len(samples)} samples')
    # In the unlikely event decorrelation math left us short (e.g. tiny counts), pad with the
    # last accepted conformer so callers always get `nsamples` structures.
    while len(samples) < nsamples:
        samples.append(coords.copy())
    return samples[:nsamples]


def build_exclusions(graph, order: int = 2) -> set:
    """Return the set of ``frozenset({i, j})`` atom-index pairs within ``order`` bonds.

    ``order=2`` excludes 1-2 (bonded) and 1-3 (shared-angle) neighbors -- the pairs held within
    hard-sphere contact by fixed bonds/angles, which must not count as overlaps.  ``graph`` is any
    object exposing ``networkx``-style ``neighbors(i)`` over integer atom indices.
    """
    import networkx as nx
    excl = set()
    lengths = dict(nx.all_pairs_shortest_path_length(graph, cutoff=order))
    for i, dists in lengths.items():
        for j, d in dists.items():
            if i < j and 1 <= d <= order:
                excl.add(frozenset((i, j)))
    return excl


def moving_set(graph, a: int, b: int) -> np.ndarray:
    """Atom indices on ``b``'s side of the ``a-b`` bond (the subtree rotated by that pivot).

    Cuts edge ``a-b`` and returns the connected component containing ``b``.  Raises if ``a-b`` is
    part of a ring (cutting it leaves ``a`` and ``b`` still connected), since a ring bond is not a
    valid independent pivot.
    """
    import networkx as nx
    g = graph.copy()
    g.remove_edge(a, b)
    comp = nx.node_connected_component(g, b)
    if a in comp:
        raise ValueError(f'bond ({a},{b}) is in a ring; not an independent rotatable bond')
    return np.array(sorted(comp), dtype=int)
