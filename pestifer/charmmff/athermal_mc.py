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
    i: int = -1
    """Anchor-side outer atom of the dihedral ``i-a-b-j`` (for the optional trans bias); ``-1`` if
    none (then this bond contributes no torsional bias)."""
    j: int = -1
    """Tip-side outer atom of the dihedral ``i-a-b-j``; ``-1`` if none."""


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


def _dihedral(p0, p1, p2, p3) -> float:
    """Signed dihedral angle (radians) of the four points p0-p1-p2-p3, in [-pi, pi]."""
    b0 = p0 - p1
    b1 = p2 - p1
    b2 = p3 - p2
    b1n = b1 / (np.linalg.norm(b1) + 1e-12)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = np.dot(v, w)
    y = np.dot(np.cross(b1n, v), w)
    return float(np.arctan2(y, x))


def _trans_penalty(phi: float) -> float:
    """Order potential ``(1 + cos phi)/2``: 0 at trans (phi=180 deg), 1 at cis.

    Biasing the acyl-tail torsions toward trans (this penalty times a tunable strength) extends
    and orders the chains -- the mechanism that carries the ensemble from fluid (Ld) toward
    liquid-ordered (Lo).  A single trans well (not the full trans+gauche rotamer set) is deliberate:
    the athermal base already samples every rotamer, and the bias only *tilts* the population toward
    extension by a tunable amount.
    """
    return 0.5 * (1.0 + np.cos(phi))


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
           max_angle: float = np.pi, seed: int = None,
           torsion_bias: float = 0.0) -> list[np.ndarray]:
    """Run dihedral-pivot MC and return ``nsamples`` decorrelated conformers.

    Starting from ``mol.coords``, propose pivots on random rotatable bonds by random angles in
    ``[-max_angle, max_angle]``; accept iff the result is overlap-free and keeps confined atoms
    inside the cylinder.  With ``torsion_bias == 0`` this is the pure athermal (T->infinity)
    criterion -- energy never enters -- and the ensemble is fluid (Ld).  With ``torsion_bias > 0``
    an extra Metropolis factor favors trans torsions (an ordering field of strength
    ``torsion_bias``), extending the chains toward the liquid-ordered (Lo) state; the caller tunes
    the strength to hit a target chain order parameter.  After ``n_equil`` proposals to melt the
    (typically rod-like) starting conformer, collect one conformer every ``n_decorr`` proposals.

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

    # Confinement is applied *monotonically inward*: a starting conformer whose tails splay wider
    # than the target cylinder (footprint > pi*R^2, common for an extended minimized lipid) begins
    # with atoms outside the wall, and a plain hard wall would reject every move and freeze the
    # start.  Instead a move is allowed while still outside if it does not push the confined atoms
    # farther out than they already are; once the whole bundle is inside R the wall is hard.  This
    # funnels an over-wide start into the cylinder and then samples within it.
    confined_idx = np.nonzero(mol.confined)[0]
    confine = np.isfinite(mol.cylinder_radius) and len(confined_idx) > 0
    if confine:
        cur_conf_max = radial_distances(coords[confined_idx], mol.axis_point,
                                        mol.axis_dir).max()

    total_proposals = n_equil + nsamples * n_decorr
    for step in range(1, total_proposals + 1):
        bond = mol.rotatable[rng.integers(nbonds)]
        theta = rng.uniform(-max_angle, max_angle)
        trial = _apply_pivot(coords, bond, theta)
        n_attempt += 1
        if confine:
            trial_conf_max = radial_distances(trial[confined_idx], mol.axis_point,
                                              mol.axis_dir).max()
            # accept iff inside the wall, or not worse than the current (still-shrinking) worst
            if trial_conf_max > mol.cylinder_radius and trial_conf_max > cur_conf_max + 1e-9:
                continue
        if _has_overlap(trial, mol.radii, bond.moving, mol.exclusions):
            continue
        if torsion_bias > 0.0 and bond.i >= 0 and bond.j >= 0:
            # Metropolis on the trans-ordering field: a move toward trans (lower penalty) is always
            # accepted; one toward cis/gauche is accepted with prob exp(-bias*dU).  bias=0 short-
            # circuits above, so the athermal path is untouched.
            du = _trans_penalty(_dihedral(trial[bond.i], trial[bond.a], trial[bond.b], trial[bond.j])) \
                - _trans_penalty(_dihedral(coords[bond.i], coords[bond.a], coords[bond.b], coords[bond.j]))
            if du > 0.0 and rng.random() >= np.exp(-torsion_bias * du):
                continue
        coords = trial
        if confine:
            cur_conf_max = trial_conf_max
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


def _element_of_mass(mass: float) -> str:
    """Coarse element symbol from atomic weight (enough to tell C/H/O/N/P/S apart in lipids)."""
    for sym, m in (('H', 1.008), ('C', 12.011), ('N', 14.007), ('O', 15.999),
                   ('P', 30.974), ('S', 32.06)):
        if abs(mass - m) < 1.0:
            return sym
    return '?'


def build_lipid_mc(coords, elements, masses, bonds, head_indices, tail_indices,
                   rmin_half, cylinder_radius, radius_scale: float = _SIGMA_FROM_RMIN,
                   min_tip_heavy: int = 2, exclusion_order: int = 3,
                   axis_dir=(0.0, 0.0, 1.0)) -> MoleculeMC:
    """Assemble a :class:`MoleculeMC` for a lipid from its topology + geometry.

    Rotatable degrees of freedom are the **acyl-tail C-C torsions**, found intrinsically: a bond
    is a tail torsion iff it is a bridge (its removal disconnects the molecule) whose tip-side
    fragment is *all carbon* -- exactly the criterion :meth:`lipid_annotate` uses to peel tails.
    This automatically excludes ester/glycerol/headgroup bonds (whose fragments carry oxygens),
    so the headgroup stays rigid and only the tails sample, per the design.  The tip-side subtree
    (hydrogens included) is what each pivot rotates.

    The confinement cylinder axis is the membrane normal (``axis_dir``, default ``+z`` -- valid
    once the molecule is head-up) through the tail bundle's lateral centroid, so radial distance
    is the in-plane footprint radius and confining the tails to ``cylinder_radius`` caps the
    footprint near the target APL.  Hard-sphere overlaps are tested only for pairs separated by
    more than ``exclusion_order`` bonds, so near-neighbor torsional (gauche) contacts -- governed
    by fixed local geometry, not excluded volume -- are not spuriously forbidden.

    Parameters
    ----------
    coords : (N,3) array_like
        Atom coordinates (the built, minimized, head-up conformer).
    elements : sequence[str] or None
        Per-atom element symbols; if ``None`` they are inferred from ``masses``.
    masses : sequence[float]
        Per-atom atomic weights (used to infer elements when ``elements`` is None).
    bonds : iterable[(int, int)]
        Bonded atom-index pairs (the full topology, hydrogens included).
    head_indices, tail_indices : sequence[int]
        Atom indices of the head reference atom(s) and the terminal tail-carbon atom(s).
    rmin_half : sequence[float]
        Per-atom CHARMM ``Rmin/2``.  Hard-sphere radii are ``rmin_half * radius_scale``.
    cylinder_radius : float
        Confinement radius (e.g. ``sqrt(APL_target/pi)``).
    radius_scale : float
        Factor on ``Rmin/2`` giving the hard-sphere radius (default scales the pair minimum down
        to the LJ collision diameter sigma).
    min_tip_heavy : int
        Minimum heavy-atom count on a bond's tip side for it to count as a rotatable torsion;
        ``2`` skips the terminal-methyl bond (a useless symmetric spin).

    Returns
    -------
    MoleculeMC
    """
    import networkx as nx

    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    if elements is None:
        elements = [_element_of_mass(m) for m in masses]
    elements = list(elements)

    full = nx.Graph()
    full.add_nodes_from(range(n))
    full.add_edges_from((int(i), int(j)) for i, j in bonds)
    heavy = full.subgraph([i for i in range(n) if elements[i] != 'H']).copy()

    rotatable: list[RotatableBond] = []
    tail_atoms: set[int] = set()
    for i, j in nx.bridges(heavy):
        g = heavy.copy()
        g.remove_edge(i, j)
        comp_i = nx.node_connected_component(g, i)
        comp_j = nx.node_connected_component(g, j)
        # the tail-tip side is the all-carbon fragment; anchor is the other side
        for tip_root, anchor_root, tip_comp in ((j, i, comp_j), (i, j, comp_i)):
            if all(elements[k] == 'C' for k in tip_comp):
                if len(tip_comp) < min_tip_heavy:
                    break
                moving = moving_set(full, anchor_root, tip_root)
                moving_ids = set(int(m) for m in moving)
                # outer dihedral atoms i-a-b-j: i on the anchor side, j on the tip side, each the
                # lowest-index carbon neighbor available (deterministic).  Used only by the optional
                # trans bias; -1 when absent.
                i_nb = _pick_neighbor(full, elements, anchor_root, exclude={tip_root})
                j_nb = _pick_neighbor(full, elements, tip_root, exclude={anchor_root},
                                      require_in=moving_ids)
                rotatable.append(RotatableBond(a=anchor_root, b=tip_root, moving=moving,
                                               i=i_nb, j=j_nb))
                tail_atoms.update(moving_ids)
                break

    confined = np.zeros(n, dtype=bool)
    if tail_atoms:
        confined[list(tail_atoms)] = True

    # Cylinder axis: the membrane normal (default +z, valid after head-up orientation), anchored
    # at the tail bundle's lateral centroid, so "radial distance" is the in-plane footprint radius
    # and confining it caps the footprint near the target APL.  Using the head->tail centroid
    # vector instead would tilt the axis toward the offset headgroup (whose ester oxygens drag the
    # centroid down) and miscenter the cylinder on the tails.
    axis_dir = np.asarray(axis_dir, dtype=float)
    axis_dir = axis_dir / np.linalg.norm(axis_dir)
    bundle = list(tail_atoms) if tail_atoms else list(range(n))
    axis_point = coords[bundle].mean(axis=0)

    radii = np.asarray(rmin_half, dtype=float) * radius_scale
    return MoleculeMC(coords=coords, radii=radii, rotatable=rotatable,
                      exclusions=build_exclusions(full, order=exclusion_order),
                      axis_point=axis_point, axis_dir=axis_dir,
                      cylinder_radius=cylinder_radius, confined=confined)


def cylinder_radius_for_apl(apl: float, inflation: float = 1.0) -> float:
    """Confinement radius for a target area-per-lipid, ``sqrt(inflation*APL/pi)``.

    A *single* lipid conformer's convex-hull footprint runs larger than its packed area-per-lipid
    (in a real bilayer neighbors interdigitate into each other's concavities), so confining the
    footprint to exactly ``APL`` would over-tighten and select for thin rods.  ``inflation`` (an
    engineering factor, calibrated so the ensemble's mean footprint lands near ``APL``) sizes the
    cylinder to ``inflation*APL``.
    """
    return float(np.sqrt(inflation * apl / np.pi))


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


def _pick_neighbor(graph, elements, node, exclude=(), require_in=None) -> int:
    """Lowest-index carbon neighbor of ``node`` (falling back to any element), for a dihedral ref.

    Excludes ``exclude`` and, if ``require_in`` is given, keeps only neighbors in that set (used to
    force the tip-side atom to lie within the rotating subtree).  Returns ``-1`` if none qualify.
    """
    exclude = set(exclude)
    cands = [k for k in graph.neighbors(node)
             if k not in exclude and (require_in is None or k in require_in)]
    if not cands:
        return -1
    cands.sort(key=lambda k: (elements[k] != 'C', k))   # carbons first, then lowest index
    return int(cands[0])


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


def tail_carbon_indices(elements, masses, bonds, n=None, min_tip_heavy: int = 2) -> list:
    """Indices of the acyl-tail carbon atoms (the aliphatic chain carbons).

    Uses the *same* bridge / all-carbon-tip criterion as :func:`build_lipid_mc`, so the tail set is
    identical to the atoms the sampler confines -- then keeps only the carbons (their hydrogens are
    dropped).  ``elements`` may be ``None``, in which case they are inferred from ``masses``.  The
    carbonyl/ester carbons and the whole headgroup are excluded automatically (their bridge tip
    fragments carry oxygens), so this returns exactly the methylene/methyl chain carbons.
    """
    import networkx as nx
    if elements is None:
        elements = [_element_of_mass(m) for m in masses]
    elements = list(elements)
    if n is None:
        n = len(elements)
    full = nx.Graph()
    full.add_nodes_from(range(n))
    full.add_edges_from((int(i), int(j)) for i, j in bonds)
    heavy = full.subgraph([i for i in range(n) if elements[i] != 'H']).copy()
    tail_atoms: set[int] = set()
    for i, j in nx.bridges(heavy):
        g = heavy.copy()
        g.remove_edge(i, j)
        comp_i = nx.node_connected_component(g, i)
        comp_j = nx.node_connected_component(g, j)
        for tip_root, anchor_root, tip_comp in ((j, i, comp_j), (i, j, comp_i)):
            if all(elements[k] == 'C' for k in tip_comp):
                if len(tip_comp) < min_tip_heavy:
                    break
                tail_atoms.update(int(m) for m in moving_set(full, anchor_root, tip_root))
                break
    return sorted(a for a in tail_atoms if elements[a] == 'C')


def chain_order_parameter(coords, elements, masses, bonds, axis_dir=(0.0, 0.0, 1.0),
                          min_tip_heavy: int = 2) -> float:
    """Tail-averaged acyl-chain order parameter of a single (head-up) conformer.

    For every C-H bond on an acyl-tail carbon (see :func:`tail_carbon_indices`), form
    ``s = 1/2 (3 cos^2(theta) - 1)`` with ``theta`` the angle between the C-H bond and the membrane
    normal ``axis_dir`` (default ``+z`` -- valid because conformers are canonicalized head-up).  The
    returned value is the tail-averaged **order** ``= -mean(s)``:

    - ``0.0``  -> isotropic / fully disordered chains (``<cos^2 theta> = 1/3``),
    - ``+0.5`` -> chains perfectly ordered along the normal (all-trans, C-H perpendicular to z).

    This is ``-S_CD`` in the usual deuterium-order-parameter convention; fluid (Ld) ensembles sit
    around 0.1-0.2 and ordered (Lo) ensembles around 0.3-0.4 (calibrate the phase targets against
    reference values).  Returns ``nan`` if the molecule has no tail C-H bonds.
    """
    import networkx as nx
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    if elements is None:
        elements = [_element_of_mass(m) for m in masses]
    elements = list(elements)
    full = nx.Graph()
    full.add_nodes_from(range(n))
    full.add_edges_from((int(i), int(j)) for i, j in bonds)
    tail_c = tail_carbon_indices(elements, None, bonds, n=n, min_tip_heavy=min_tip_heavy)
    axis = np.asarray(axis_dir, dtype=float)
    axis = axis / np.linalg.norm(axis)
    svals = []
    for c in tail_c:
        for h in full.neighbors(c):
            if elements[h] == 'H':
                v = coords[h] - coords[c]
                nv = np.linalg.norm(v)
                if nv == 0.0:
                    continue
                cval = float(np.dot(v / nv, axis))
                svals.append(0.5 * (3.0 * cval * cval - 1.0))
    if not svals:
        return float('nan')
    return float(-np.mean(svals))


def ensemble_chain_order(conformers, elements, masses, bonds, axis_dir=(0.0, 0.0, 1.0),
                         min_tip_heavy: int = 2) -> float:
    """Mean :func:`chain_order_parameter` over a list of conformer coordinate arrays.

    ``conformers`` is an iterable of ``(N,3)`` coordinate arrays that share the one ``bonds`` /
    ``elements`` topology (the usual case: many samples of one molecule).  ``nan`` conformers (no
    tail C-H) are skipped; returns ``nan`` if none contribute.
    """
    vals = [chain_order_parameter(c, elements, masses, bonds, axis_dir, min_tip_heavy)
            for c in conformers]
    vals = [v for v in vals if v == v]  # drop nan
    return float(np.mean(vals)) if vals else float('nan')
