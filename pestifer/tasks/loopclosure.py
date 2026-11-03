# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Closing modeled interior loops onto their downstream anchors: CCD closure, clearing any ring
piercing the closure made, and bonding the junction (``connect``).

A mixin, because two tasks run it.  :class:`~pestifer.tasks.ligate.LigateTask` did, after a minimize;
since 2026-10-30 :class:`~pestifer.tasks.psfgen.PsfgenTask` does, right after building, whenever a
``ligate`` task follows -- so the first minimize sees every loop closed, clash-filtered and bonded,
instead of ``guesscoord``'s placement.  That placement can run straight through resolved protein:
example 9's modeled gp41 F 29-43A loop was built through chain A 69-76 (0.72 A at closest), and the
minimize that pulled them apart left cis bonds in both, which a later closure -- re-seeding only phi
and psi -- carried to the final structure.  ``ligate`` then finds its loops already closed and stands
aside; its ``ccd`` settings are the ones psfgen used.
"""
import logging
import shutil

from ..core.artifacts import *
from ..molecule.molecule import Molecule
from ..scripters import GenericScripter, PsfgenScripter

logger = logging.getLogger(__name__)

_RECLOSE_ATTEMPTS = 3
"""How many times :meth:`LoopClosureMixin.reclose_piercing_loops` re-closes a loop that still pierces
a ring after rotation, each time from a fresh seed, before leaving it to terminate's gate."""
_PEPTIDE_BOND = 1.33
_JUNCTION_BOND_SLACK = 0.5
"""How far (A) a CCD-closed loop's C may sit from the ideal peptide bond to its anchor's N."""
_JUNCTION_MIN_OMEGA = 150.0
"""The least trans (|omega|, deg) a CCD-closed loop's junction may be."""


class LoopClosureMixin:
    """CCD closure, piercing repair and ``connect`` for a task holding ``base_molecule``."""

    _ccd_context = None
    """What :meth:`close_loops_ccd` closed and how, kept so a loop can be re-closed later."""

    #: prefix of this closure's log lines, so a psfgen-time closure does not read as ligate's
    _closure_tag = 'ligate'

    def with_obstacles(self, extra_env=None):
        """``extra_env`` plus the heavy atoms of ``source.obstacles``, for a modeler's clash
        environment.  Every loop and tail closure goes through here, so a part modeled with
        another molecule excluded still keeps clear of where that molecule will be."""
        import numpy as np
        mol = getattr(self, 'base_molecule', None)
        obs = mol.obstacle_coords() if mol is not None else None
        parts = [np.asarray(x, dtype=float) for x in (extra_env, obs) if x is not None and len(x)]
        return np.vstack(parts) if parts else None

    def close_and_connect_loops(self, specs: dict) -> int:
        """Close every interior loop by CCD (``specs['ccd']``), clear ring piercings the closure
        made (unless ``check_piercings`` is false), bond the junctions, and check stereochemistry
        at the new junctions.  ``specs`` is a ``ligate`` task's.  Returns 0 on success."""
        logger.debug('Closing loops onto their downstream anchors by cyclic coordinate descent')
        result = self.close_loops_ccd(specs.get('ccd', {}) or {})
        if result != 0:
            return result
        return self.finish_loop_closure(specs)

    def finish_loop_closure(self, specs: dict) -> int:
        """After the loops are placed (by CCD or by steering): clear ring piercings, bond the
        junctions, check stereochemistry at them.  Marks ``base_molecule.loops_closed_by``."""
        # BEFORE `connect`, and that ordering is the whole thing.  While the loop is still
        # unpatched its C-terminus is free, so the branch carrying a threading bond is a small
        # pendant hanging off a backbone bridge -- rotatable.  `connect` joins the loop at BOTH
        # ends, and from then on any hinge inside it would swing the entire rest of the chain,
        # so `pendant_axes` skips every axis (`max_pendant` is 180 atoms) and reports "no
        # rotatable loop hinge found".  Measured on example 11 seed 27021974: with the call
        # placed after `connect` the resolver detected `HSD D-53 <- ALA D-50`, could not clear
        # it, warned, and the build ran ~30 more minutes before `terminate` refused to package
        # it.  A resolver that cannot succeed is worse than none: it reports the problem in a
        # voice that sounds like it handled it.
        if specs.get('check_piercings', True):
            self.resolve_loop_piercings()
        logger.debug('Connecting loop C-termini to their partner N-termini')
        result = self.connect(specs.get('connect', {}))
        if result != 0:
            return result
        self.base_molecule.loops_closed_by = self.basename
        # The closures' junctions exist only from here; checked now, a cis bond or D residue made
        # by the closure is blamed on the closing task, not on the minimize that follows it.
        self.track_stereochemistry()
        return 0
    def resolve_loop_piercings(self):
        """Detect and clear any ring piercing caused by a model-built loop, by rotating the
        offending loop sub-branch out of the ring.

        Closing a loop is an **unphysical displacement**: CCD drives backbone dihedrals to meet
        a target, and the `steer` method pulls termini together with a colvars bias.  Both put
        in the energy needed to cross a van der Waals barrier, which free dynamics at a kT of
        about 0.6 cannot -- so both can thread a bond through a ring, exactly as a minimizer or
        a grid placement can.  The clash-filtered ensemble cannot see one: a threading has no
        atomic clash, which is the whole difficulty.

        Nothing looked before this.  Nine shipped examples (7-12, 14, 15, 18) run `ligate` and
        then `md` with no check in between, and until the terminate gate landed nothing checked
        the final structure either.  Scanning the packages those nine produced in the 3.22.1
        sweep found zero piercings, so this is a mechanism being closed rather than an observed
        failure being fixed -- but a piercing is topological and survives every later stage, so
        a clean final structure is real evidence only about those runs.

        Mirrors :meth:`PsfgenTask.resolve_tail_piercings`: same engine
        (:func:`~pestifer.psfutil.ring_resolve.resolve_model_built_piercings`), same
        best-effort, non-fatal contract, toggled with ``ligate.check_piercings``.  A piercing
        no rotation clears is reported and left to a downstream ``ring_check``.

        **Runs before ``connect``, while the loop is still a pendant.**  An unpatched modeled
        loop has a free C-terminus, so the branch carrying a threading bond hangs off a
        backbone bridge and can be swung out.  Once ``connect`` patches the loop to its
        downstream anchor the loop is joined at both ends, and a hinge inside it would carry
        the whole rest of the chain -- :meth:`RingChecker.pendant_axes` skips any branch over
        ``max_pendant`` (180) atoms, so every axis is rejected and nothing can be resolved.

        Serials are read from the current state PDB rather than captured earlier, because the
        closure itself rewrites coordinates; ``loop_atoms_from_pdb`` looks the loop up by
        segname and resid, so it does not depend on serial numbering being stable.
        """
        from ..psfutil.loop_ccd import loop_atoms_from_pdb
        from ..psfutil.ring_resolve import resolve_model_built_piercings

        gaps = self.base_molecule.protein_loop_gaps()
        if not gaps:
            return
        state: StateArtifacts = self.get_current_artifact('state')
        loop_serials = set()
        for gap in gaps:
            try:
                _o, _c, ser = loop_atoms_from_pdb(state.pdb.name, gap['loop_resids'],
                                                  segname=gap['segname'])
            except Exception as e:
                logger.debug(f'{self._closure_tag}: could not locate loop {gap["segname"]} '
                             f'{gap["loop_resids"][0]}-{gap["loop_resids"][-1]} in the '
                             f'post-connect structure ({type(e).__name__}: {e})')
                continue
            loop_serials.update(ser)
        if not loop_serials:
            return
        xsc = state.xsc.name if getattr(state, 'xsc', None) else None
        self.next_basename('loop-piercing-check')
        fixed_pdb, unresolved = resolve_model_built_piercings(
            state.psf.name, state.pdb.name, xsc, self.basename, sorted(loop_serials),
            segtypes=['protein', 'glycan'], mover_kind='loop')
        if fixed_pdb is not None:
            shutil.copy(fixed_pdb, f'{self.basename}.pdb')
            self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf,
                               xsc=getattr(state, 'xsc', None)),
                          key='state', artifact_type=StateArtifacts)
            logger.info(f'{self._closure_tag}: resolved loop ring-piercing(s) -> {self.basename}.pdb')
        # A rotation inside a pendant loop swings its free end too, off the anchor the closure put
        # it on, and connect would then bond across the gap (ex11 D:52-53 at 6.8 A, closed cis).
        # Only a CCD closure put the end on an ideal junction; a steered one never sat there.
        opened = self.open_junctions(gaps) if self._ccd_context is not None else {}
        if (unresolved or opened) and self._ccd_context is not None:
            unresolved, opened = self.reclose_piercing_loops(unresolved, loop_serials, opened)
        for i, (d, w) in opened.items():
            g = gaps[i]
            logger.warning(f'{self._closure_tag}: loop {g["segname"]}:{g["loop_resids"][0]}-'
                           f'{g["loop_resids"][-1]} no longer meets its anchor {g["segname"]}:'
                           f'{g["c_anchor_resid"]} (C-N {d:.2f} A, omega {w:.0f} deg) after piercing '
                           f'resolution; connect bonds it as it lies, and the junction may close cis')
        if unresolved:
            how = ('rotation or by re-closing the loop' if self._ccd_context is not None
                   else 'rotation')
            logger.warning(f'{len(unresolved)} loop ring-piercing(s) could not be cleared by '
                           f'{how}; terminate will refuse to package the result unless a later '
                           f'ring_check clears it.  Add a ring_check task after the next '
                           f'minimize, or rebuild with a different --seed:')
            for p in unresolved:
                ee, er = p['piercee'], p['piercer']
                logger.warning(f'  ring of {ee.get("resname", "?")} {ee["segname"]}-{ee["resid"]} '
                               f'pierced by a bond in loop {er["segname"]}-{er["resid"]}')

    def open_junctions(self, gaps):
        """Loops in the current state whose free end no longer meets its downstream anchor:
        ``{gap index: (C-N distance, omega)}`` for each junction more than
        :data:`_JUNCTION_BOND_SLACK` off a peptide bond or less trans than :data:`_JUNCTION_MIN_OMEGA`.

        CCD puts the last loop residue's CA and C on an ideal trans peptide upstream of the anchor,
        so a closed junction sits near 1.33 A and 180 deg; one moved after closure does not."""
        from ..psfutil.loop_ccd import backbone_from_pdb, dihedral_deg
        import numpy as np
        out = {}
        try:
            pdb = self.get_current_artifact('state').pdb.name
            bbs = {}
            for i, g in enumerate(gaps):
                seg = g['segname']
                if seg not in bbs:
                    bbs[seg] = backbone_from_pdb(pdb, segname=seg)
                last, anc = bbs[seg].get(g['loop_resids'][-1]), bbs[seg].get(g['c_anchor_resid'])
                if not last or not anc:
                    continue
                d = float(np.linalg.norm(last['C'] - anc['N']))
                w = float(dihedral_deg(last['CA'], last['C'], anc['N'], anc['CA']))
                if abs(d - _PEPTIDE_BOND) > _JUNCTION_BOND_SLACK or abs(w) < _JUNCTION_MIN_OMEGA:
                    out[i] = (d, w)
        except Exception as e:
            logger.debug(f'{self._closure_tag}: junction check skipped ({type(e).__name__}: {e})')
        return out

    def reclose_piercing_loops(self, unresolved, loop_serials, opened=None):
        """Re-close, from a fresh seed, each loop that still pierces a ring after rotation --
        this time with every other loop present, as placed.

        :meth:`close_loops_ccd` closes all loops at once, each against the resolved structure
        only, so no loop can see another.  Two copies of one long loop closed that way can
        interpenetrate, and a threading is not a clash, so nothing in the closure notices.
        Example 14 replica 1 in the 3.29.8 sweep: the symmetry copies of the 37-residue loop
        657-693 closed into each other, a bond of A SER 678 through the ring of E TYR 675.  No
        rotation could clear it, and terminate refused the structure 1.5 hours of MD later.

        What is measured is that a fresh seed clears it: on that structure, re-closing loop A
        cleared the piercing on the first attempt, and so did five other seeds both with and
        without the other loops present.  The other loops are put in the environment anyway,
        since a pose that overlaps them is one the least-clash selection should reject; that it
        matters has not been shown.  The loop is re-closed if the piercing bond lies in it;
        failing that, the loop holding the pierced ring.  Each attempt is followed by the same
        rotation pass and check.

        A loop whose junction a rotation opened (``opened``, from :meth:`open_junctions`) is
        re-closed the same way, and every rotation pass after a re-closure is checked for that
        too, since it can open another.

        Returns ``(unresolved, opened)``: the piercings and open junctions still left after
        :data:`_RECLOSE_ATTEMPTS` attempts.
        """
        from ..psfutil.loop_ccd import loop_atoms_from_pdb, close_one_loop, heavy_coords_of_serials
        from ..psfutil.ring_resolve import resolve_model_built_piercings
        from ..util.coord import pdb_replace_coords
        import numpy as np

        ctx = self._ccd_context
        gaps = ctx['gaps']
        all_loop_serials = ctx['all_loop_serials']

        def gap_holding(res):
            for i, g in enumerate(gaps):
                if g['segname'] == res.get('segname') and res.get('resid') in g['loop_resids']:
                    return i
            return None

        opened = dict(opened or {})
        reclosed = []
        for attempt in range(1, _RECLOSE_ATTEMPTS + 1):
            targets = set(opened)
            for p in unresolved:
                i = gap_holding(p['piercer'])
                if i is None:
                    i = gap_holding(p['piercee'])
                if i is not None:
                    targets.add(i)
            if not targets:
                break
            state: StateArtifacts = self.get_current_artifact('state')
            src_pdb = state.pdb.name
            new_coords = {}
            for i in sorted(targets):
                g = gaps[i]
                tag = f"{g['segname']}:{g['loop_resids'][0]}-{g['loop_resids'][-1]}"
                _o, _c, own = loop_atoms_from_pdb(src_pdb, g['loop_resids'], segname=g['segname'])
                others = sorted(set(all_loop_serials) - set(own))
                seed = ctx['seed'] + attempt * len(gaps) + i
                res = close_one_loop(src_pdb=src_pdb, segname=g['segname'],
                                     loop_resids=g['loop_resids'],
                                     n_anchor_resid=g['n_anchor_resid'],
                                     c_anchor_resid=g['c_anchor_resid'],
                                     all_loop_serials=all_loop_serials, seed=seed,
                                     extra_env=self.with_obstacles(heavy_coords_of_serials(src_pdb, others)),
                                     **ctx['kwargs'])
                for k, serial in enumerate(res['serials']):
                    new_coords[serial] = res['closed'][k]
                reclosed.append({'loop': tag, 'attempt': attempt, 'seed': seed})
                logger.info(f'{self._closure_tag}/ccd: re-closed loop {tag} against the other loops as placed '
                            f'(attempt {attempt} of {_RECLOSE_ATTEMPTS}, seed {seed})')
            self.next_basename(f'reclose-{attempt}')
            out_pdb = f'{self.basename}.pdb'
            serial_list = sorted(new_coords)
            pdb_replace_coords(src_pdb, out_pdb, np.array([new_coords[s] for s in serial_list]),
                               {s: k for k, s in enumerate(serial_list)})
            self.pdb_to_coor(out_pdb)
            self.register(dict(psf=state.psf, pdb=PDBFileArtifact(self.basename),
                               coor=NAMDCoorFileArtifact(self.basename)),
                          key='state', artifact_type=StateArtifacts)
            xsc = state.xsc.name if getattr(state, 'xsc', None) else None
            fixed_pdb, unresolved = resolve_model_built_piercings(
                state.psf.name, out_pdb, xsc, self.basename, sorted(loop_serials),
                segtypes=['protein', 'glycan'], mover_kind='loop')
            if fixed_pdb is not None:
                self.next_basename(f'reclose-{attempt}-rotated')
                shutil.copy(fixed_pdb, f'{self.basename}.pdb')
                self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf,
                                   xsc=getattr(state, 'xsc', None)),
                              key='state', artifact_type=StateArtifacts)
            opened = self.open_junctions(gaps)
            if not unresolved and not opened:
                logger.info(f'{self._closure_tag}/ccd: re-closing cleared every loop ring-piercing '
                            f'and left every loop on its anchor ({attempt} attempt(s))')
                break
        self.record_outcome(loop_reclosures=reclosed, loop_piercings_unresolved=len(unresolved),
                            loop_junctions_open=len(opened))
        return unresolved, opened

    def close_loops_ccd(self, ccd_specs: dict) -> int:
        """
        Close each modeled interior loop onto its downstream anchor by cyclic coordinate
        descent (offline, deterministic). For each gap the loop's C-terminal backbone is
        walked onto an ideal peptide-bond target built from the anchor, distributing the
        closure across the loop's own phi/psi dihedrals (whole residues, sidechains included).

        These gaps are floppy, solvent-exposed surface loops (that is *why* they are
        unresolved): the goal is not to recover a unique native conformation but to produce a
        *physically plausible, closed, clash-free* starting structure for MD. Each loop is
        seeded from realistic Ramachandran basins (not a straight chain) and, for each gap, a
        small **ensemble** of independent seeds is closed and the candidate that introduces the
        fewest steric clashes with the rest of the structure is kept -- preserving the one
        virtue of the old straight-chain+steer approach (it added no clashes) while giving a
        far more defensible backbone. The closed coordinates replace the loop atoms in a new
        state PDB; ``connect`` then patches the LINK bond and rebuilds the junction, and a light
        minimization (if enabled) relaxes residual strain.
        """
        import os
        import sys
        import numpy as np
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
        from scipy.spatial.distance import cdist
        from ..psfutil.loop_ccd import loop_atoms_from_pdb, close_one_loop
        from ..util.coord import pdb_replace_coords
        from ..util.progress import PestiferProgress

        self.next_basename('ccd')
        state: StateArtifacts = self.get_current_artifact('state')
        src_pdb = state.pdb.name
        gaps = self.base_molecule.protein_loop_gaps()
        if not gaps:
            logger.warning(f'{self._closure_tag} method=ccd: has_protein_loops is set but no interior gaps '
                           'were enumerated; nothing to close')
            return 0
        # `ligate.ccd.seed` if set, else the build's namd.seed (BaseTask.derived_seed).  It used to
        # default to a constant `--seed` never reached, so every replica closed its loops from the
        # same seed -- and since the closure's input is a non-reproducible minimize, replicas
        # diverged anyway, without the seed meaning anything.
        seed = self.derived_seed(ccd_specs.get('seed'))
        self.record_outcome(ccd_seed=seed)
        max_iters = int(ccd_specs.get('max_iters', 2000))
        tol = float(ccd_specs.get('tol', 0.1))
        ensemble = max(1, int(ccd_specs.get('ensemble', 10)))
        refine_iters = int(ccd_specs.get('refine', 250))
        guard = bool(ccd_specs.get('guard', False))
        on_clash = str(ccd_specs.get('on_clash', 'warn')).casefold()
        kwargs = dict(ensemble=ensemble, refine=refine_iters, guard=guard, max_iters=max_iters,
                      tol=tol)

        # Serials of EVERY modeled loop copy: while unclosed they hold throwaway guesscoord
        # positions, so all are excluded from every loop's closure environment.
        all_loop_serials = set()
        for gap in gaps:
            _o, _c, gser = loop_atoms_from_pdb(src_pdb, gap['loop_resids'], segname=gap['segname'])
            all_loop_serials.update(gser)
        all_loop_serials = sorted(all_loop_serials)
        self._ccd_context = dict(gaps=gaps, all_loop_serials=all_loop_serials, seed=seed,
                                 kwargs=kwargs)

        def _args(i, extra_env=None):
            g = gaps[i]
            return dict(src_pdb=src_pdb, segname=g['segname'], loop_resids=g['loop_resids'],
                        n_anchor_resid=g['n_anchor_resid'], c_anchor_resid=g['c_anchor_resid'],
                        all_loop_serials=all_loop_serials, seed=seed + i, ensemble=ensemble,
                        refine=refine_iters, guard=guard, max_iters=max_iters, tol=tol,
                        extra_env=self.with_obstacles(extra_env))

        # Loops that do not interfere can be closed at once (steering closes all loops
        # simultaneously too). Close every loop CONCURRENTLY against the resolved structure only;
        # a downstream minimize resolves whatever mutual (non-threaded) clashes result. The one
        # thing minimize cannot undo is one loop's backbone threaded through another's, so a
        # post-pass detects loops whose backbones interpenetrate and RE-CLOSES those sequentially
        # (each avoiding the ones before it).
        #
        # **Per-loop seed = base + index, so this is deterministic GIVEN ITS INPUT -- which is
        # not the same as a reproducible build, and the difference matters.** `ligate` is
        # preceded by a minimize, and parallel NAMD is not bitwise reproducible: measured
        # 2026-10-21, the same `.namd` on byte-identical input at +p24 twice gave two different
        # `.coor` files. So CCD gets a different starting geometry every run and closes the
        # loops differently, while the seed controls nothing a user can see.
        #
        # That is how example 11 came to pierce a ring in sweep replica 3 and not in 1 or 2: until
        # 3.29.6 `--seed` never reached `ccd.seed`, so all three closed from the identical seed and
        # still diverged.  The seed now derives from `--seed`, but the input still varies run to
        # run, so even a fixed seed does not reproduce a closure.
        #
        # Do not read the word "deterministic" here as "rerunning reproduces the closure". It
        # does not, and a rerun that comes back clean is evidence of nothing.
        n = len(gaps)
        results = [None] * n
        nworkers = max(1, min(n, (os.cpu_count() or 2)))

        # Determinate progress bar with a time-remaining estimate. The total work is a known
        # count -- n loops x `ensemble` candidate closures -- so `PestiferProgress(max_value=...)`
        # gives a real ETA rather than a spinner. Each worker pushes one tick per candidate onto
        # `q`; the driver drains it and, when a loop finishes, trues-up progress to that loop's full
        # ensemble share (a loop can break early after fewer than `ensemble` candidates). Only shown
        # on an interactive terminal so batch logs stay clean.
        total = n * max(1, ensemble)
        prog = {'done': 0, 'loops': 0}
        show = sys.stdout.isatty()
        bar = None
        if show:
            bar = PestiferProgress(name=f'{self._closure_tag}/ccd', color='fuchsia', max_value=total)
            bar.register_update_function(lambda: min(1.0, prog['done'] / total) if total else 1.0)

        def _drain(q):
            while True:
                try:
                    q.get_nowait()
                except Exception:
                    break
                prog['done'] += 1

        if n == 1 or nworkers == 1:
            # serial path: no mid-call draining is possible, so progress advances per loop
            for i in range(n):
                results[i] = close_one_loop(**_args(i))
                prog['loops'] += 1
                prog['done'] = prog['loops'] * ensemble
                if bar:
                    bar.go()
        else:
            with mp.Manager() as mgr:
                q = mgr.Queue()
                with ProcessPoolExecutor(max_workers=nworkers) as ex:
                    futs = {ex.submit(close_one_loop, progress_queue=q, **_args(i)): i
                            for i in range(n)}
                    pending = set(futs)
                    while pending:
                        completed, pending = wait(pending, timeout=0.25,
                                                  return_when=FIRST_COMPLETED)
                        _drain(q)
                        for fut in completed:
                            results[futs[fut]] = fut.result()
                            prog['loops'] += 1
                        # true-up: a finished loop always counts its full ensemble share
                        prog['done'] = max(prog['done'], prog['loops'] * ensemble)
                        if bar:
                            bar.go()
                    _drain(q)
        if bar:
            prog['done'] = total
            bar.finish()
        logger.debug(f'{self._closure_tag}/ccd: closed {n} loop(s) in parallel across {nworkers} worker(s)')

        # Optional threading repair (OFF by default). Closing all loops at once lets copies
        # converging at an assembly axis overlap; on the 4zmj trimer those overlaps were ordinary
        # clashes that the downstream minimize pushes apart, so the parallelism is kept. Set
        # ccd.thread_ca > 0 to re-close, sequentially, any pair whose backbones come closer than
        # that Ca-Ca distance.
        #
        # The interlinked case this comment once called "not observed" was observed: example 14
        # replica 1 in the 3.29.8 sweep, one copy of loop 657-693 threaded through a TYR ring of
        # the other. A Ca-Ca threshold is the wrong trigger for that (the piercing bond was a side
        # chain), so `reclose_piercing_loops` re-closes on the piercing itself, after rotation.
        THREAD_CA = float(ccd_specs.get('thread_ca', 0.0))
        if THREAD_CA > 0:
            adj = [[] for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    if cdist(results[i]['ca'], results[j]['ca']).min() < THREAD_CA:
                        adj[i].append(j); adj[j].append(i)
            seen = [False] * n
            for i in range(n):
                if seen[i]:
                    continue
                comp, stack = [], [i]
                while stack:
                    u = stack.pop()
                    if seen[u]:
                        continue
                    seen[u] = True; comp.append(u)
                    stack.extend(v for v in adj[u] if not seen[v])
                if len(comp) > 1:
                    comp.sort()
                    logger.info(f'{self._closure_tag}/ccd: {len(comp)} loops interpenetrate after the parallel '
                                f'closure; re-closing them sequentially (ccd.thread_ca={THREAD_CA})')
                    placed = [results[comp[0]]['heavy']]   # keep the first; re-close the rest
                    for idx in comp[1:]:
                        results[idx] = close_one_loop(**_args(idx, extra_env=np.vstack(placed)))
                        placed.append(results[idx]['heavy'])

        # collect coordinates and emit per-loop steric diagnostics
        new_coords = {}
        broken = []
        rotation_report = []      # (tag, segname, rotation rows) per loop, for the provenance file
        for i, g in enumerate(gaps):
            res = results[i]
            rep = res['rep']
            tag = f"{g['segname']}:{g['loop_resids'][0]}-{g['loop_resids'][-1]} (len {len(g['loop_resids'])})"
            rotation_report.append((tag, g['segname'], res.get('rotations', [])))
            ndeep = rep['n_deep'] + rep['n_env_deep']
            if rep['min_ca'] < 3.0:     # a crossed backbone -- minimization cannot undo it
                broken.append((tag, rep))
                logger.warning(
                    f'{self._closure_tag}/ccd: loop {tag} is THREADED (min non-adjacent CA {rep["min_ca"]:.2f} A, '
                    f'worst overlap {rep["worst"]:.2f} A) -- a crossed backbone that minimization '
                    f'cannot undo. Increase ccd.ensemble/refine, or provide an externally-modeled '
                    f'loop (see docs/design/loop-modeling.md).')
            elif ndeep:
                logger.info(
                    f'{self._closure_tag}/ccd: loop {tag} closed with {ndeep} deep + '
                    f'{rep["n_soft"] + rep["n_env_soft"]} soft overlaps (worst {rep["worst"]:.2f} A, '
                    f'min CA {rep["min_ca"]:.2f} A) -- not threaded; a downstream minimize relaxes these.')
            elif rep['n_soft'] or rep['n_env_soft']:
                logger.info(f'{self._closure_tag}/ccd: loop {tag} closed with soft contacts only '
                            f'({rep["n_soft"] + rep["n_env_soft"]}; worst {rep["worst"]:.2f} A).')
            for k, serial in enumerate(res['serials']):
                new_coords[serial] = res['closed'][k]

        if broken and on_clash == 'error':
            logger.error(f'{self._closure_tag}/ccd: {len(broken)} loop(s) threaded; aborting (ccd.on_clash=error)')
            return 1

        self.write_rotation_report(rotation_report, kind='loop closure')

        out_pdb = f'{self.basename}.pdb'
        serial_list = sorted(new_coords)
        coords_arr = np.array([new_coords[s] for s in serial_list])
        row_of_serial = {s: i for i, s in enumerate(serial_list)}
        pdb_replace_coords(src_pdb, out_pdb, coords_arr, row_of_serial)
        self.pdb_to_coor(out_pdb)
        self.register(dict(psf=state.psf,
                           pdb=PDBFileArtifact(self.basename, pytestable=True),
                           coor=NAMDCoorFileArtifact(self.basename)),
                      key='state', artifact_type=StateArtifacts)
        return 0
    
    def connect(self, connect_specs):
        """
        Connect the loop termini using the specified ``LINK`` patch.
        """
        logger.debug(f'Connect specs: {connect_specs} (unused)')
        self.write_connect_patches()
        result = self.connect_gaps()
        return result

    def write_connect_patches(self):
        """
        Write the connect patches to a data file.
        This method generates a data file that contains the connection patches for the loop termini.
        It uses the base molecule to write the connection patches and updates the state variables with the new data file.
        The data file is named based on the current basename, which is generated by the ``next_basename`` method.
        """
        self.next_basename('gap_patches')
        mol: Molecule = self.base_molecule
        datafile = f'{self.basename}.inp'
        writer: GenericScripter = self.get_scripter('data')
        writer.newfile(datafile)
        mol.write_connect_patches(writer)
        writer.writefile()
        self.register(self.basename, key='connect_patches_input', artifact_type=InputFileArtifact)

    def connect_gaps(self):
        """
        Connect the gaps in the loop termini.
        This method uses the psfgen scripter to create a script that connects the gaps in the loop termini
        using the specified patch. It generates a new PSF file and PDB file based on the current state of the base molecule
        and the connection patches defined in the data file. The script is then executed, and if successful,
        the resulting PSF and PDB files are saved in the current state.
        
        Returns
        -------
        int
            The result of the psfgen script execution. A return value of 0 indicates success, while any other value indicates failure.
        """
        self.next_basename('heal')
        pg: PsfgenScripter = self.get_scripter('psfgen')
        pg.newscript(self.basename)
        RM: ResourceManager = self.resource_manager
        CC: CHARMMFFContent = RM.charmmff_content
        CC.copy_charmmfile_local('pestifer.top')
        charmm_topology_files: FileArtifactList = self.get_current_artifact_data('charmmff_topfiles')
        charmm_topology_files.append(CharmmffTopFileArtifact('pestifer', ext='top'))
        pg.addline(f'topology pestifer.top')
        patchfile: Path = self.get_current_artifact_path('connect_patches_input')
        state: StateArtifacts = self.get_current_artifact('state')
        pg.load_project(state.psf.name, state.pdb.name)
        pg.addline(f'source {patchfile.name}')
        pg.writescript(self.basename, guesscoord=True, regenerate=True)
        self.register(self.basename, key='connect_gaps_tcl', artifact_type=PsfgenInputScriptArtifact)
        result = pg.runscript()
        if result == 0:
            self.pdb_to_coor(f'{self.basename}.pdb')
            self.register(dict(
                psf=PSFFileArtifact(self.basename, pytestable=True), 
                pdb=PDBFileArtifact(self.basename, pytestable=True), 
                coor=NAMDCoorFileArtifact(self.basename)), key='state', artifact_type=StateArtifacts)
            self.register(self.basename, key='connect_gaps_log', artifact_type=LogFileArtifact)
        return result
