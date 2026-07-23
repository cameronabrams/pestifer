# Author: Cameron F. Abrams <cfa22@drexel.edu>
"""
Physics-based modeling of model-built protein *terminal tails* (free N-/C-terminal runs).

Companion to :mod:`pestifer.psfutil.loop_ccd`. An interior gap has two resolved anchors and
must be *closed* onto the downstream one (CCD); a terminal tail has only one anchor and a free
end, so there is nothing to close -- the job is to replace ``guesscoord``'s arbitrary extended
arm with a *realistic, clash-free* backbone that still projects cleanly away from the fold.
This regularizes terminal-tail modeling to the same standard the interior loop closer gives
(``docs/design/loop-modeling.md``): a Ramachandran-seeded backbone plus a clash-filtered
ensemble, differing only in that closure is replaced by a rigid **junction placement**.

Per tail, per ensemble member:

1. **Sample** each tail residue's (phi, psi) from the coil Ramachandran basins and apply them
   to the built full-atom tail (sidechains carried) -- a locally realistic *shape*, in an
   arbitrary pose.
2. **Place the junction** -- Kabsch-superpose the three backbone atoms of the tail's anchored
   residue onto an ideal trans-peptide target built off the resolved anchor, carrying the whole
   tail rigidly. The free end is unconstrained (correctly), so no closure is needed.
3. **Score** the pose's steric overlap with the frozen environment (:func:`loop_clash_report`)
   and keep the least-clashing member.

Everything operates on numpy arrays; deterministic for a fixed seed. A downstream ``minimize``
relaxes residual soft overlaps, exactly as for interior loops.
"""
import numpy as np

from ..util.coord import kabsch, rotate_points_about_axis
from .loop_ccd import (
    backbone_from_pdb, loop_atoms_from_pdb, build_loop_problem, apply_backbone_dihedrals,
    sample_backbone_dihedrals, place_atom_nerf, anchor_closure_target,
    heavy_env_coords_from_pdb, loop_clash_report, _clash_score,
    extract_backbone, loop_rotation_report,
)

# Standard trans-peptide internal coordinates for building a junction target.
_PEP = dict(C_N=1.33, N_CA=1.45, CA_C=1.52,
            CA_C_N=116.6, C_N_CA=121.7, N_CA_C=111.0)


def downstream_anchor_target(anchor_N, anchor_CA, anchor_C, psi_deg=150.0, phi_deg=-120.0):
    """
    Ideal target positions for a C-terminal tail's *first* residue backbone atoms ``(N, CA, C)``
    such that that residue forms a good trans peptide bond *downstream* of the resolved anchor
    (``anchor:C -> tail_first:N``).

    This is the mirror of :func:`pestifer.psfutil.loop_ccd.anchor_closure_target` (which targets
    the residue *upstream* of an anchor N, for interior-loop closure). Built from the fixed anchor
    N/CA/C by NeRF placement with standard peptide geometry; ``psi_deg`` (anchor N-CA-C-N torsion)
    and ``phi_deg`` (the first tail residue's C-prev-N-CA-C torsion) are the two free junction
    parameters -- a generic extended value is fine because the tail's free end carries no further
    constraint and a downstream ``minimize`` relaxes the exact junction. Returns ``(N, CA, C)``.
    """
    N = place_atom_nerf(anchor_N, anchor_CA, anchor_C, _PEP['C_N'], _PEP['CA_C_N'], psi_deg)
    CA = place_atom_nerf(anchor_CA, anchor_C, N, _PEP['N_CA'], _PEP['C_N_CA'], 180.0)
    C = place_atom_nerf(anchor_C, N, CA, _PEP['CA_C'], _PEP['N_CA_C'], phi_deg)
    return np.array([N, CA, C])


def _place_junction(shaped, junction_rows, target):
    """Rigid-body-superpose the anchored-residue backbone atoms onto ``target`` and carry the
    whole tail; returns the placed coordinates."""
    R, t, _ = kabsch(shaped[junction_rows], target)
    return shaped @ R.T + t


def _refine_tail(shaped, prob, order, tail_resids, junction_rows, target, env, rng,
                 n_iters, perturb_deg=25.0, T0=8.0, Tmin=0.4):
    """
    Iteratively declash a placed tail by perturb-and-replace simulated annealing -- the free-tail
    analogue of :func:`pestifer.psfutil.loop_ccd.refine_declash_ccd` (which perturbs-and-recloses).

    Each step perturbs one backbone dihedral, re-places the junction by rigid superposition (a
    tail has no closure target, so re-anchoring replaces re-closure), and scores steric overlap
    with the frozen environment. Moves lowering the clash score are accepted; worse moves pass a
    Metropolis test under a geometrically annealed temperature. An occasional large swing about a
    leading (long-lever) bond rotates the tail bulk (e.g. away from a crowded axis); local moves
    polish. Deterministic for a fixed ``rng``.
    """
    bonds, masks = prob['bonds'], prob['moving_masks']
    nb = len(bonds)

    def evaluate(shaped_pose):
        placed = _place_junction(shaped_pose, junction_rows, target)
        rep = loop_clash_report(order, placed, tail_resids, env_coords=env)
        return _clash_score(rep), rep, placed

    cur = np.array(shaped, dtype=float)
    cur_s, cur_rep, cur_placed = evaluate(cur)
    best, best_s, best_rep, best_placed = np.array(cur), cur_s, cur_rep, cur_placed
    if best_s == 0:
        return best_placed, best_rep
    early = max(1, nb // 3)                       # leading third = longest levers
    for it in range(n_iters):
        T = T0 * (Tmin / T0) ** (it / max(1, n_iters - 1))
        if rng.random() < 0.35:
            b = int(rng.integers(early)); delta = float(rng.normal(0.0, 90.0))
        else:
            b = int(rng.integers(nb)); delta = float(rng.normal(0.0, perturb_deg))
        a, c = bonds[b]
        mask = np.asarray(masks[b], dtype=bool)
        trial = np.array(cur)
        trial[mask] = rotate_points_about_axis(trial[mask], trial[a], trial[c] - trial[a], delta)
        s, rep, placed = evaluate(trial)
        if s <= cur_s or rng.random() < np.exp(-(s - cur_s) / max(T, 1e-6)):
            cur, cur_s = trial, s
            if s < best_s:
                best_s, best_rep, best_placed = s, rep, placed
                if best_s == 0:
                    break
    return best_placed, best_rep


def model_one_tail(src_pdb, segname, tail_resids, anchor_resid, end,
                   all_modeled_serials, seed, ensemble=10, refine=200, extra_env=None):
    """
    Model one built terminal tail against the resolved structure -- the standalone unit the
    psfgen declash step runs per tail (parallel to
    :func:`pestifer.psfutil.loop_ccd.close_one_loop`).

    Reads the tail and its environment from ``src_pdb``: the environment is every heavy atom NOT
    in ``all_modeled_serials`` (all model-built loops/tails hold throwaway coords) and NOT in the
    tail's own anchor (the bonded junction residue). Samples each residue's Ramachandran (phi,
    psi), places the anchored junction by rigid superposition, and keeps the least-clashing of
    ``ensemble`` seeds. Deterministic for a fixed ``seed``.

    Parameters
    ----------
    src_pdb : str
    segname : str
        Actual (image-mapped) segid of this tail copy.
    tail_resids : sequence[int]
        The model-built tail residues, N->C.
    anchor_resid : int
        The adjoining resolved residue (fixed junction).
    end : {'N', 'C'}
        Which terminus the tail extends. 'N' -> the tail is upstream of the anchor (its last
        residue bonds anchor:N); 'C' -> downstream (its first residue's N bonds anchor:C).
    all_modeled_serials : iterable[int]
        PDB serials of every model-built loop/tail atom (all excluded from the environment).
    seed : int
    ensemble : int
    refine : int
        Iterative-declash refinement iterations per candidate (perturb-and-replace simulated
        annealing); 0 = ensemble only. A downstream minimize relaxes any residual soft overlaps.
    extra_env : (K, 3) array, optional
        Extra heavy-atom coordinates to treat as fixed environment -- e.g. tails already modeled
        this pass, so converging tails (a trimer's C-termini near the assembly axis) declash
        against one another instead of interpenetrating.

    Returns
    -------
    dict : ``segname``, ``tail`` (resids), ``serials``, ``order``, ``modeled`` (M,3),
    ``rep`` (:func:`loop_clash_report`).
    """
    rng = np.random.default_rng(seed)
    bb = backbone_from_pdb(src_pdb, segname=segname)
    order, coords, serials = loop_atoms_from_pdb(src_pdb, tail_resids, segname=segname)
    prob = build_loop_problem(order, coords, tail_resids)
    row = prob['row']

    # environment: everything heavy except all modeled loops/tails and this tail's anchor,
    # plus any already-modeled tails passed as extra_env.
    _ao, _ac, anch_ser = loop_atoms_from_pdb(src_pdb, [anchor_resid], segname=segname)
    env = heavy_env_coords_from_pdb(src_pdb, exclude_serials=list(all_modeled_serials) + list(anch_ser))
    if extra_env is not None and len(extra_env):
        extra_env = np.asarray(extra_env, dtype=float)
        env = np.vstack([env, extra_env]) if len(env) else extra_env

    # which of the anchored residue's backbone atoms the junction target constrains. The target
    # itself is rebuilt per candidate below with a *sampled* junction torsion, because the
    # direction the tail emanates from the anchor is a real DOF: with it fixed, the base residue
    # is pinned and no internal-torsion move can pull it out of a clash (the failure mode
    # observed on the crowded gp41 C-termini).
    aN, aCA, aC = bb[anchor_resid]['N'], bb[anchor_resid]['CA'], bb[anchor_resid]['C']
    if end == 'N':
        jr = tail_resids[-1]                      # anchored (upstream) residue bonds anchor:N
        junction_rows = [row[(jr, 'CA')], row[(jr, 'C')], row[(jr, 'O')]]
        make_target = lambda tor: anchor_closure_target(aN, aCA, aC, phi_deg=tor)   # (CA, C, O)
    else:
        jr = tail_resids[0]                        # anchored (downstream) residue's N bonds anchor:C
        junction_rows = [row[(jr, 'N')], row[(jr, 'CA')], row[(jr, 'C')]]
        make_target = lambda tor: downstream_anchor_target(aN, aCA, aC, psi_deg=tor)  # (N, CA, C)

    # A synthetic prev-C reference so the first residue's phi is well-defined; the free end's
    # absolute torsion is arbitrary (the pose is fixed by the junction Kabsch below).
    n0 = bb[tail_resids[0]]['N']
    ca0 = bb[tail_resids[0]]['CA']
    prev_C = n0 + (n0 - ca0)

    best = None
    for _cand in range(max(1, ensemble)):
        phipsi = sample_backbone_dihedrals(rng, len(tail_resids))
        shaped = apply_backbone_dihedrals(prob['coords'], prob, tail_resids, phipsi,
                                          prev_C=prev_C, next_N=n0)   # next_N unused (last psi free)
        target = make_target(float(rng.uniform(-180.0, 180.0)))     # sampled emanation direction
        if refine > 0:
            placed, rep = _refine_tail(shaped, prob, order, tail_resids, junction_rows, target,
                                       env, rng, n_iters=refine)
        else:
            placed = _place_junction(shaped, junction_rows, target)
            rep = loop_clash_report(order, placed, tail_resids, env_coords=env)
        key = (rep['n_deep'] + rep['n_env_deep'], rep['n_soft'] + rep['n_env_soft'])
        if best is None or key < best[0]:
            best = (key, placed, rep)
        if key[0] == 0 and key[1] == 0:
            break

    _key, modeled, rep = best
    heavy_rows = [i for i, (_rr, nm) in enumerate(order) if not nm.startswith('H')]
    # provenance: the internal backbone rotations from guesscoord (bb) to the modeled handoff.
    # The anchored end supplies the junction reference; the free end has no phi/psi partner.
    prev_C, next_N = (None, aN) if end == 'N' else (aC, None)
    rotations = loop_rotation_report(bb, extract_backbone(order, modeled), tail_resids,
                                     prev_C=prev_C, next_N=next_N)
    return dict(segname=segname, tail=list(tail_resids), serials=serials, order=order,
                modeled=modeled, rep=rep, heavy=modeled[heavy_rows], rotations=rotations)
