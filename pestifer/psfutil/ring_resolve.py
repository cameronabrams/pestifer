# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Shared ring-piercing *resolution* by glycan-pendant rotation.

A grafted (or otherwise model-built) glycan bond can thread through the hole of a rigid ring --
a proline side-chain ring, another glycan's pyranose ring, or an aromatic protein side-chain ring
-- without any atoms overlapping.  Such a threading is *topological*: it has a low/zero
atomic-clash count, so a clash-count declash cannot see it, yet it survives minimization and
RATTLE-fails at the first dynamics step.

The fix is to rotate the offending glycan sub-branch about an upstream rotatable bond (a molecular
bridge) so the piercing bond is pulled back out of the ring.  This module holds that resolution
logic so it can be used both by the up-front :class:`~pestifer.tasks.ringcheck.RingCheckTask`
(which also handles aromatic-ring-side rotations and lipid deletion) and at glycan *graft time*
inside the psfgen task (:meth:`PsfgenTask.declash`), where a threaded glycan is rotated out
before the input structure is written -- ideally making a downstream ``ring_check`` unnecessary.

The rotation moves coordinates only (never topology): candidates are scored in memory and the
pose that un-threads the ring with the fewest new heavy-atom clashes is written.
"""
import logging

import numpy as np

from .psfring import ring_check, RingChecker
from ..util.coord import rotate_points_about_axis, pdb_replace_coords
from ..util.util import cell_from_xsc
from ..util.contact_gate import DEFAULT_MIN_CONTACT as HARD_OVERLAP

logger = logging.getLogger(__name__)

#: Heavy-atom separation (A) a rotation may not CREATE.  Tied to the pre-flight steric
#: gate's threshold: a pose inside it is refused before dynamics anyway, so accepting it
#: only moves the failure later and loses the chance to try another hinge.  See
#: `pestifer.util.contact_gate` for the measurements behind the number.

# pendant rotations tried, in order, to pull a glycan piercing bond out of a rigid ring
_PENDANT_DEGREES = [60, 120, 180, 240, 300, 90, 270]


def _fmt(side: dict) -> str:
    return f'{side.get("resname", "?")} {side["segname"]}-{side["resid"]}'


def _is_glycan_pierce(p: dict) -> bool:
    """True if piercing ``p`` is a glycan bond threading a protein/glycan ring -- the case the
    glycan-pendant rotation can clear (by moving the glycan, never the protein)."""
    piercee, piercer = p['piercee'], p['piercer']
    return (piercer['segtype'] == 'glycan' and bool(piercer.get('bond_serials'))
            and piercee['segtype'] in ('protein', 'glycan'))


def try_glycan_pendant(checker: RingChecker, working_pdb: str, base: np.ndarray, box,
                       piercer: dict, target: tuple, out_pdb: str, piercee_label: str = '',
                       degrees=None, mover_serials=None, mover_kind: str = 'glycan'):
    """Rotate the sub-branch carrying the piercing bond about an upstream rotatable bond.

    Hinge axes come from :meth:`RingChecker.pendant_axes` (smallest branch first); each is swept
    over ``degrees``.  Among the rotations that clear the ring, the one with the fewest new
    heavy-atom clashes wins and only that single pose is written to ``out_pdb``.  The rotation is
    done in memory (no per-candidate PDB).

    Despite the name the rotation is chemistry-agnostic (it rigidly rotates a graph branch); pass
    ``mover_serials`` to restrict the eligible hinges to those whose moving branch stays *within*
    a given atom set -- e.g. a model-built terminal tail, so only the tail swings and the resolved
    chain is never moved. ``mover_kind`` names the branch in the log message.

    Returns ``out_pdb`` on success or ``None`` if no rotation clears the piercing.
    """
    if degrees is None:
        degrees = _PENDANT_DEGREES

    _NEAREST = 10          # pendant_axes' own default; named so the escalation can say what it widened FROM

    def _axes(max_axes):
        ax = checker.pendant_axes(piercer['bond_serials'], piercer['segname'],
                                  max_axes=max_axes)
        if mover_serials is None:
            return ax
        want = set(int(s) for s in mover_serials)
        serial_of_row = {r: s for s, r in checker._row_of_serial.items()}
        return [a for a in ax if all(serial_of_row[int(r)] in want for r in a['prows'])]

    def _search(axes):
        """Best pose over these hinges, or None. (clash, coords, i_ser, j_ser, deg)"""
        best = None
        scratch = base.copy()
        for ax in axes:
            prows = ax['prows']
            pivot = base[ax['j_row']]
            axis = base[ax['j_row']] - base[ax['i_row']]
            branch = base[prows]
            for deg in degrees:
                scratch[prows] = rotate_points_about_axis(branch, pivot, axis, deg)
                if checker.check_coords(scratch, box, [target]):
                    continue  # still pierced
                # ELIGIBILITY BEFORE RANKING.  "Fewest clashes wins" is not the same as "few
                # clashes": when every candidate is bad it keeps the least bad and hands it on.
                # Measured on an example-17 patchB, the accepted poses scored 3, 3, 3, 4, 0, 2,
                # 21 and 7 new contacts, and the structure went to dynamics with two heavy atoms
                # 0.334 A apart.  A pose that puts a heavy-atom pair inside the steric-gate
                # threshold is not a resolution -- the gate downstream refuses it anyway, so
                # accepting it only moves the failure later and loses the chance to try another
                # hinge.  Filter first, then rank the survivors by the softer count.
                if checker.clash_count(scratch, prows, cutoff=HARD_OVERLAP, box=box):
                    continue
                clash = checker.clash_count(scratch, prows, box=box)
                if best is None or clash < best[0]:
                    best = (clash, scratch.copy(), ax['i_ser'], ax['j_ser'], deg)
                if clash == 0:
                    break
            scratch[prows] = branch  # restore before trying the next axis
            if best is not None and best[0] == 0:
                break
        return best

    axes = _axes(_NEAREST)
    if not axes:
        logger.debug(f'  no rotatable {mover_kind} hinge found for {_fmt(piercer)}')
        return None
    best = _search(axes)
    if best is None:
        # ESCALATE.  `pendant_axes` returns hinges smallest-moving-branch-first and keeps only
        # the first few, which optimizes for least collateral motion -- and that is exactly the
        # wrong half of the list for un-threading a tail.  Measured on job 26469605's stuck
        # piercings (2026-10-13): production tried 10 hinges with moving branches of 27-52
        # atoms, while the hinge that clears the piercing moves ~130 of the lipid's ~134 --
        # near the headgroup, swinging the whole tail.  All five stuck piercings across both
        # stages clear from the full hinge set, with ZERO new clashes, using these same angles.
        # The cap, not the move set and not the angle resolution, is what sent 4-5 cholesterols
        # per build to deletion.
        #
        # Widening only on failure costs nothing on the piercings that already resolve, and is
        # paid only where the alternative is deleting a molecule.
        tried = {(a['i_ser'], a['j_ser']) for a in axes}
        wider = [a for a in _axes(None) if (a['i_ser'], a['j_ser']) not in tried]
        if wider:
            logger.debug(f'  {_fmt(piercer)}: no hinge among the nearest {len(axes)} cleared it; '
                         f'widening to {len(wider)} more')
            best = _search(wider)
    if best is None:
        return None
    clash, coords, i_ser, j_ser, deg = best
    pdb_replace_coords(working_pdb, out_pdb, coords, checker._row_of_serial)
    logger.info(f'  {piercee_label or _fmt({"segname": target[0], "resid": target[1]})}: '
                f'rotating the {_fmt(piercer)} {mover_kind} branch {deg} deg about bond '
                f'{i_ser}-{j_ser} clears the piercing ({clash} new heavy-atom clash(es))')
    return out_pdb


#: Sterol long-axis spin angles.  Coarse on purpose: the acceptance test is "does it un-thread",
#: so a fine sweep only costs time.  Mirrors the staging of _PENDANT_DEGREES -- the large,
#: least-biased rotations first, then the quarter turns.
_SPIN_DEGREES = [180, 90, 270, 120, 240, 60, 300, 150, 210, 30, 330]


def try_sterol_spin(checker: RingChecker, working_pdb: str, base: np.ndarray, box,
                    piercee: dict, target: tuple, out_pdb: str, piercer_label: str = '',
                    degrees=None):
    """Rigidly spin the pierced molecule about its own long axis to un-thread it.

    The partner move to :func:`try_glycan_pendant`. That one bends the piercing tail; this one
    turns the speared ring. Either can clear a lipid-tail-through-sterol-ring piercing, and they
    fail in different situations: a tail hinged near its tip has little leverage, while a sterol
    wedged between neighbors cannot turn.

    The axis is the molecule's first principal axis -- for a sterol, the ring plate plus its
    iso-octyl tail -- taken through its centroid, so the rotation preserves the molecule's
    position and its membrane-normal alignment and only turns it in place. That matters: a sterol
    displaced laterally would trade one piercing for a clash.

    Honest about what this can and cannot do: a fused ring system is roughly coaxial with its own
    long axis, so spinning it moves the ring ATOMS much more than it moves the HOLE. It clears a
    piercing whose bond passes near the rim, and will not clear one threaded through the centre.
    The caller gets ``None`` for the latter and falls through to :func:`try_piercee_translate`.

    An earlier version of this docstring added that "a sterol displaced laterally would trade one
    piercing for a clash", and used that to justify having no translation move at all.  Measured
    2026-10-22 on example 17's stuck CHL1 L1-197: of 1022 rigid translations, 636 un-thread the
    ring and **25 do so with no hard overlap**, the best at 1.5 A for 5 soft contacts.  The claim
    was wrong, and it was load-bearing -- it is why a sterol skewered through two rings had no
    remedy but deletion.

    Returns ``out_pdb`` on success or ``None`` if no spin clears the piercing.
    """
    if degrees is None:
        degrees = _SPIN_DEGREES
    seg, resid = piercee['segname'], str(piercee['resid'])
    rows = np.array([checker._row_of_serial[a.serial] for a in checker._atoms
                     if a.segname == seg and str(a.resid.resid) == resid], dtype=int)
    if rows.size < 3:
        logger.debug(f'  {_fmt(piercee)}: too few atoms to define a spin axis')
        return None
    pts = base[rows]
    centroid = pts.mean(axis=0)
    # first principal axis = the molecule's long axis
    _u, _s, vt = np.linalg.svd(pts - centroid, full_matrices=False)
    axis = vt[0]
    best = None
    scratch = base.copy()
    for deg in degrees:
        scratch[rows] = rotate_points_about_axis(pts, centroid, axis, deg)
        if checker.check_coords(scratch, box, [target]):
            continue  # still pierced
        # same eligibility rule as `_search`: a spin that buries the sterol in a neighbour is
        # not a resolution, and the steric gate refuses it downstream anyway
        if checker.clash_count(scratch, rows, cutoff=HARD_OVERLAP, box=box):
            continue
        clash = checker.clash_count(scratch, rows, box=box)
        if best is None or clash < best[0]:
            best = (clash, scratch.copy(), deg)
        if clash == 0:
            break
    if best is None:
        return None
    clash, coords, deg = best
    pdb_replace_coords(working_pdb, out_pdb, coords, checker._row_of_serial)
    logger.info(f'  {_fmt(piercee)}: spinning it {deg} deg about its long axis clears the '
                f'piercing by {piercer_label or "a bond"} ({clash} new heavy-atom clash(es))')
    return out_pdb


def _piercer_within(p: dict, mover_serials: set) -> bool:
    """True if piercing ``p`` is caused by a bond whose atoms all lie in ``mover_serials`` -- i.e.
    a model-built branch (loop/tail) we are free to rotate."""
    bs = p['piercer'].get('bond_serials')
    return bool(bs) and set(int(s) for s in bs) <= mover_serials


#: Rigid-translation search for :func:`try_piercee_translate`, smallest first.  A piercing is
#: cleared by the least displacement that works, so these are tried in order and the search stops
#: at the first magnitude that yields any eligible pose.
_TRANSLATE_MAGNITUDES = (0.5, 1.0, 1.5, 2.0, 2.5, 3.0)

#: Number of directions sampled on the unit sphere.  A deterministic Fibonacci lattice, not an
#: RNG draw: a resolver that gives a different answer on a rerun cannot be debugged, and this
#: repository has already paid for that lesson once in `ligate`'s loop closure.
_TRANSLATE_DIRECTIONS = 98


def _unit_sphere(n: int) -> np.ndarray:
    """``n`` roughly-uniform directions on the unit sphere (Fibonacci lattice, deterministic)."""
    k = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * k
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


def try_piercee_translate(checker: RingChecker, working_pdb: str, base: np.ndarray, box,
                          piercee: dict, target: tuple, out_pdb: str, piercer_label: str = '',
                          magnitudes=None):
    """Slide the whole pierced molecule rigidly until the threading bond is outside its ring.

    The last remedy, after :func:`try_glycan_pendant` (bend the piercing tail) and
    :func:`try_sterol_spin` (turn the speared ring).  Both of those are ROTATIONS, and a
    rotation cannot free a molecule that is skewered at two points: example 17's CHL1 L1-197 was
    threaded through its A-ring (6-membered, the 3-OH end) by one acyl chain and through its
    D-ring (5-membered, the tail end) by another, each crossing within 0.2 A of the ring centre.
    A sterol's long axis runs near both ring centres, so no spin moves them, and bending either
    chain leaves the other one holding it.  Every one of 903 and 917 candidate tail rotations
    left the ring threaded.  A translation moves both rings away from both bonds at once: 1.5 A
    cleared it for 5 soft contacts.

    **Only a molecule may be translated, never a residue inside one.**  Moving one residue of a
    covalent chain would stretch the bonds either side of it into nonsense, so this refuses
    unless the residue's atoms have no bonds leaving the residue.  In a bilayer a lipid is its
    own molecule and the guard passes; for a protein or glycan ring it refuses and the caller
    reports the piercing instead, which is correct -- a glycan threaded onto a protein is a
    report about the stage that built it, not something to shove sideways.

    Searched smallest-displacement-first over a deterministic direction set, with the same
    eligibility rule as the rotations: a pose that puts a heavy-atom pair inside the steric gate
    is not a resolution.  Among the eligible poses at the first magnitude that yields any, the
    fewest new soft contacts wins.

    Returns ``out_pdb`` on success or ``None``.
    """
    if magnitudes is None:
        magnitudes = _TRANSLATE_MAGNITUDES
    seg, resid = piercee['segname'], str(piercee['resid'])
    serials = [a.serial for a in checker._atoms
               if a.segname == seg and str(a.resid.resid) == resid]
    if not serials:
        logger.debug(f'  {_fmt(piercee)}: no atoms found to translate')
        return None

    # the molecule guard -- see the docstring.  A residue bonded to anything outside itself is
    # part of a larger covalent structure and must not be moved on its own.
    own = set(serials)
    G = checker.topol.G
    outside = {nb for s in serials if s in G for nb in G.neighbors(s)} - own
    if outside:
        logger.debug(f'  {_fmt(piercee)}: not a free molecule ({len(outside)} bond(s) leave the '
                     f'residue), so a rigid translation would break covalent geometry; refusing')
        return None

    rows = np.array([checker._row_of_serial[s] for s in serials], dtype=int)
    pts = base[rows]
    dirs = _unit_sphere(_TRANSLATE_DIRECTIONS)
    scratch = base.copy()
    for mag in magnitudes:
        best = None
        for d in dirs:
            scratch[rows] = pts + mag * d
            if checker.check_coords(scratch, box, [target]):
                continue  # still pierced
            if checker.clash_count(scratch, rows, cutoff=HARD_OVERLAP, box=box):
                continue  # buries it in a neighbour; the steric gate refuses this downstream
            clash = checker.clash_count(scratch, rows, box=box)
            if best is None or clash < best[0]:
                best = (clash, scratch.copy(), d.copy())
        if best is not None:
            clash, coords, d = best
            pdb_replace_coords(working_pdb, out_pdb, coords, checker._row_of_serial)
            logger.info(f'  {_fmt(piercee)}: sliding it {mag:.1f} A '
                        f'({d[0]:+.2f} {d[1]:+.2f} {d[2]:+.2f}) clears the piercing by '
                        f'{piercer_label or "its partner"} ({clash} new heavy-atom clash(es))')
            return out_pdb
    scratch[rows] = pts
    logger.debug(f'  {_fmt(piercee)}: no rigid translation up to {magnitudes[-1]} A clears it')
    return None


def resolve_model_built_piercings(psf: str, pdb: str, xsc: str, out_prefix: str, mover_serials,
                                  *, cutoff: float = 3.5, segtypes=('protein', 'glycan'),
                                  max_ring_size: int = 7, mover_kind: str = 'tail'):
    """Detect and clear ring piercings caused by a *model-built* backbone branch (a terminal tail
    or loop), by rotating only that branch out of the ring.

    Analogous to :func:`resolve_glycan_piercings`, but the resolvable set is defined by atom
    membership (the piercing bond lies wholly within ``mover_serials``) rather than by segtype, and
    the hinge is restricted so only those atoms move (the resolved structure is never touched). A
    model-built tail has a free end, so swinging its distal portion about a backbone bridge is
    unconstrained -- the natural fit. Piercings *not* caused by ``mover_serials`` are ignored here
    (they are another pass's concern), so only mover-caused piercings appear in ``unresolved``.

    Returns ``(fixed_pdb, unresolved)`` as :func:`resolve_glycan_piercings` does.
    """
    mover_serials = set(int(s) for s in mover_serials)
    piercings = ring_check(psf, pdb, xsc, cutoff=cutoff, segtypes=list(segtypes),
                           max_ring_size=max_ring_size)
    resolvable = [p for p in piercings if _piercer_within(p, mover_serials)]
    if not resolvable:
        return None, []
    checker = RingChecker(psf, cutoff=cutoff, segtypes=list(segtypes), max_ring_size=max_ring_size)
    box = cell_from_xsc(xsc)[0] if xsc else None
    working_pdb = pdb
    unresolved = []
    for p in resolvable:
        piercee = p['piercee']
        target = (piercee['segname'], piercee['resid'])
        base = checker.load_coords(working_pdb)
        out_pdb = f'{out_prefix}-{mover_kind}-{piercee["segname"]}-{piercee["resid"]}.pdb'
        cand = try_glycan_pendant(checker, working_pdb, base, box, p['piercer'], target,
                                  out_pdb, _fmt(piercee), mover_serials=mover_serials,
                                  mover_kind=mover_kind)
        if cand is None:
            unresolved.append(p)
        else:
            working_pdb = cand
    fixed_pdb = working_pdb if working_pdb != pdb else None
    return fixed_pdb, unresolved


def resolve_glycan_piercings(psf: str, pdb: str, xsc: str, out_prefix: str, *,
                             cutoff: float = 3.5, segtypes=('protein', 'glycan'),
                             max_ring_size: int = 7):
    """Detect ring piercings caused by glycan bonds and clear each by rotating the glycan pendant.

    Only piercings whose *piercer* is a glycan bond (threading a protein or glycan ring) are
    resolvable this way -- the glycan is moved, the protein is not.  Piercings this path cannot
    address (e.g. a ring speared by a lipid bond) are returned as unresolved for the caller to
    handle (or ignore).

    Parameters
    ----------
    psf, pdb, xsc : str
        The structure to check/fix; ``xsc`` may be ``None`` (non-periodic / vacuum).
    out_prefix : str
        Prefix for the (single) written fix PDB(s), e.g. ``<basename>``.
    cutoff, segtypes, max_ring_size
        Forwarded to the piercing scan / :class:`RingChecker`.

    Returns
    -------
    (fixed_pdb, unresolved) : (str | None, list[dict])
        ``fixed_pdb`` is the path of the final resolved coordinates, or ``None`` if nothing was
        changed (no glycan piercings, or none clearable).  ``unresolved`` lists the piercings that
        remain (glycan piercings no rotation cleared, plus any non-glycan-piercer piercings).
    """
    piercings = ring_check(psf, pdb, xsc, cutoff=cutoff, segtypes=list(segtypes),
                           max_ring_size=max_ring_size)
    if not piercings:
        return None, []
    resolvable = [p for p in piercings if _is_glycan_pierce(p)]
    unresolved = [p for p in piercings if not _is_glycan_pierce(p)]
    if not resolvable:
        return None, unresolved

    checker = RingChecker(psf, cutoff=cutoff, segtypes=list(segtypes), max_ring_size=max_ring_size)
    box = cell_from_xsc(xsc)[0] if xsc else None
    working_pdb = pdb
    for p in resolvable:
        piercee = p['piercee']
        target = (piercee['segname'], piercee['resid'])
        base = checker.load_coords(working_pdb)
        out_pdb = f'{out_prefix}-glycan-{piercee["segname"]}-{piercee["resid"]}.pdb'
        cand = try_glycan_pendant(checker, working_pdb, base, box, p['piercer'], target,
                                  out_pdb, _fmt(piercee))
        if cand is None:
            unresolved.append(p)
        else:
            working_pdb = cand
    fixed_pdb = working_pdb if working_pdb != pdb else None
    return fixed_pdb, unresolved
