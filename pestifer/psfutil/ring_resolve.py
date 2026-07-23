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

logger = logging.getLogger(__name__)

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
    axes = checker.pendant_axes(piercer['bond_serials'], piercer['segname'])
    if mover_serials is not None:
        mover_serials = set(int(s) for s in mover_serials)
        serial_of_row = {r: s for s, r in checker._row_of_serial.items()}
        axes = [ax for ax in axes
                if all(serial_of_row[int(r)] in mover_serials for r in ax['prows'])]
    if not axes:
        logger.debug(f'  no rotatable {mover_kind} hinge found for {_fmt(piercer)}')
        return None
    best = None  # (clash, coords, i_ser, j_ser, deg)
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
            clash = checker.clash_count(scratch, prows)
            if best is None or clash < best[0]:
                best = (clash, scratch.copy(), ax['i_ser'], ax['j_ser'], deg)
            if clash == 0:
                break
        scratch[prows] = branch  # restore before trying the next axis
        if best is not None and best[0] == 0:
            break
    if best is None:
        return None
    clash, coords, i_ser, j_ser, deg = best
    pdb_replace_coords(working_pdb, out_pdb, coords, checker._row_of_serial)
    logger.info(f'  {piercee_label or _fmt({"segname": target[0], "resid": target[1]})}: '
                f'rotating the {_fmt(piercer)} {mover_kind} branch {deg} deg about bond '
                f'{i_ser}-{j_ser} clears the piercing ({clash} new heavy-atom clash(es))')
    return out_pdb


def _piercer_within(p: dict, mover_serials: set) -> bool:
    """True if piercing ``p`` is caused by a bond whose atoms all lie in ``mover_serials`` -- i.e.
    a model-built branch (loop/tail) we are free to rotate."""
    bs = p['piercer'].get('bond_serials')
    return bool(bs) and set(int(s) for s in bs) <= mover_serials


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
