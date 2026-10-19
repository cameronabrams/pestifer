# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`RingCheckTask` class for checking for pierced rings in a molecular structure.
This class is a descendant of the :class:`BaseTask <pestifer.tasks.basetask.BaseTask>` class and is used to check for pierced rings in a molecular structure.
It identifies configurations where a ring is pierced by a bond and resolves each by whichever motion is available:

- an **aromatic** protein side-chain ring (His/Phe/Tyr/Trp) pivots on its own dihedral, so it is swung off the piercing bond by rotating the side chain (chi2, then chi1);
- a **rigid** ring that cannot itself move -- a **proline** side-chain ring (fused to the backbone) or a **glycan** ring -- when speared by a **glycan** bond is cleared by rotating the *glycan pendant* so the piercing bond is pulled out of the ring;
- a **glycan** ring speared by a **protein side chain** is cleared by rotating the *glycan branch
  carrying the pierced ring* off the side chain (the glycan is the model-built partner, so it moves
  first); failing that, the piercing **side chain** itself is rotated about its chi bonds;
- a pierced **lipid** ring deletes a segment, chosen by ``delete``: ``piercee``, ``piercer``, ``both`` (the piercee plus the piercer when it too is a lipid, for inter-threaded lipids), or ``none`` (report and stop);
- anything that no available rotation clears stops the build with the residue named.

Which segtypes are checked is set by ``segtypes`` (there is no default, so it must be given explicitly).
Usage is described in the :ref:`config_ref tasks ring_check` documentation.
"""
import glob
import logging
import os
import shutil

import numpy as np

from .basetask import BaseTask
from ..scripters import PsfgenScripter
from ..core.artifacts import *
from ..core.errors import PestiferBuildError
from ..psfutil.psfring import ring_check, RingChecker
from ..psfutil.ring_resolve import HARD_OVERLAP, try_glycan_pendant, try_sterol_spin
from ..util.util import cell_from_xsc

logger=logging.getLogger(__name__)

#: Rotation passes before giving up.  Clearing a piercing can create another, so one
#: pass is not enough; a pass that fails to reduce the count will not do better the
#: next time either, so the loop also stops on no net progress.
_MAX_RESOLVE_CYCLES = 5

# aromatic protein residues whose side-chain ring pivots freely on chi2/chi1 (proline is
# NOT here: its ring is fused into the backbone, so a chi rotation cannot swing it clear)
_AROMATIC = {'HIS', 'HSD', 'HSE', 'HSP', 'PHE', 'TYR', 'TRP'}
# Backbone atom names.  Everything else in a residue is its side chain, and restricting the
# hinge search to that set makes the bond graph pick out exactly the residue's chi bonds -- no
# per-residue chi table, correct for every residue, and it naturally finds nothing for GLY (no
# side chain) or PRO (side chain fused to the backbone, so no chi bond is a graph bridge).
_BACKBONE = {'N', 'HN', 'H', 'HT1', 'HT2', 'HT3', 'CA', 'HA', 'HA1', 'HA2',
             'C', 'O', 'OT1', 'OT2', 'OXT'}
# side-chain dihedral rotations tried, in order, to swing an aromatic ring off the bond
_CHI2_DEGREES = [180, 120, 240, 90, 270, 60, 300]
_CHI1_DEGREES = [120, 240, 180, 90, 270]


class RingCheckTask(BaseTask):
    """
    A class for checking for pierced rings
    """

    _yaml_header = 'ring_check'
    """
    YAML header for RingCheckTask objects.
    This header is used to declare RingCheckTask objects in YAML task lists.
    """

    @staticmethod
    def _fmt(side):
        return f'{side.get("resname", "?")} {side["segname"]}-{side["resid"]}'

    def _report(self, piercings, level):
        for p in piercings:
            level(f'  ring of {self._fmt(p["piercee"])} pierced by a bond in '
                  f'{self._fmt(p["piercer"])} ({p["piercer"]["segtype"]})')

    @staticmethod
    def _rotatable(p):
        """True if some rotation can be tried on this piercing: an aromatic side-chain ring
        (rotate the ring) or any rigid protein/glycan ring speared by a glycan bond (rotate
        the glycan pendant)."""
        piercee, piercer = p['piercee'], p['piercer']
        aromatic = piercee['segtype'] == 'protein' and piercee.get('resname') in _AROMATIC
        pendant = (piercee['segtype'] in ('protein', 'glycan')
                   and piercer['segtype'] == 'glycan' and piercer.get('bond_serials'))
        # a glycan ring speared by a protein side chain: either partner can move
        sidechain = (piercee['segtype'] == 'glycan' and piercer['segtype'] == 'protein'
                     and (piercee.get('ring_serials') or piercer.get('bond_serials')))
        # A lipid ring speared by a lipid bond -- in practice a sterol ring threaded by an acyl
        # tail, which is the overwhelming majority of membrane piercings.  These used to be
        # excluded, so they never reached the rotation engine at all and went straight to
        # deletion: an ex17 build deleted 7-11 lipids per run without attempting a single
        # rotation.  At one or two piercings deletion is tolerable; at ten it changes the
        # composition of the membrane being built.
        lipid = (piercee['segtype'] == 'lipid' and piercer['segtype'] == 'lipid'
                 and (piercer.get('bond_serials') or piercee.get('ring_serials')))
        # A lipid ring speared by a PROTEIN bond -- a sterol threaded by a side chain of an
        # embedded protein.  This was the one combination with no entry here, so it fell through
        # to the lipid deletion path and `delete: piercee` removed the CHOLESTEROL while leaving
        # the side chain doing the threading exactly where it was: the sterol depletion d89cc518
        # eliminated, re-entering through a different door, and not even a fix.
        #
        # Measured, not anticipated: `CHL1 159 ring <- LEU 702 bond` on the replica of the
        # 2026-10-16 example-17 arm that died post-embed.  It is invisible to that build because
        # ex17 runs `segtypes: [lipid, glycan]`, so nothing had ever reached this branch.
        lipid_protein = (piercee['segtype'] == 'lipid' and piercer['segtype'] == 'protein'
                         and (piercer.get('bond_serials') or piercee.get('ring_serials')))
        return bool(aromatic or pendant or sidechain or lipid or lipid_protein)

    def do(self):
        state: StateArtifacts = self.get_current_artifact('state')
        cutoff: float = self.specs.get('cutoff', 3.5)
        segtypes: list = self.specs.get('segtypes', []) or []
        delete_these: str = self.specs.get('delete', 'piercee')
        max_ring_size: int = self.specs.get('max_ring_size', 7)
        resolve: bool = bool(self.specs.get('resolve', True))
        if not segtypes:
            # Was a warning that returned 0 -- a `ring_check` task that checked nothing and
            # reported success, which is the silent-no-op shape this repository keeps a list of.
            # A threaded ring is a topologically wrong molecule that no later dynamics can undo,
            # so a check configured to look at nothing is a configuration error.  Reachable only
            # by writing `segtypes: []` on purpose, since an omitted key now takes the default.
            raise PestiferBuildError(
                'ring_check: `segtypes` is empty, so this task would check nothing and report '
                'success. Omit the key to get the default [lipid, glycan, protein], or name the '
                'segtypes you want explicitly.')

        def run_check(st):
            xsc = st.xsc.name if st.xsc else None
            if xsc is None:
                logger.debug('No XSC in pipeline state — ring_check runs in non-periodic (vacuum) mode')
            return ring_check(st.psf.name, st.pdb.name, xsc, cutoff=cutoff,
                              segtypes=segtypes, max_ring_size=max_ring_size)

        npiercings = run_check(state)
        n_initial = len(npiercings)
        n_rotations = 0

        # VERIFY-ONLY.  `resolve: false` is the terminator of the minimize-after-rotation loop
        # spliced by MakeMembraneSystemTask._guard_piercings: by this point the structure has
        # had its rotations relaxed by a minimize, and anything still threaded is something two
        # rounds of rotate-and-relax could not clear.  Rotating here would hand dynamics a
        # structure whose coordinates were last changed by a rigid move nothing minimized --
        # which is the exact failure this pass exists to prevent -- so it reports and stops.
        if not resolve:
            self.record_outcome(piercings_initial=n_initial, rotations=0,
                                piercings_remaining=n_initial, resolve=False)
            if not npiercings:
                self.result = 0
                return self.result
            self._report(npiercings, logger.error)
            named = '; '.join(f'{self._fmt(p["piercee"])} pierced by {self._fmt(p["piercer"])}'
                              for p in npiercings)
            raise PestiferBuildError(
                f'{len(npiercings)} ring piercing(s) survived the rotate-and-minimize loop: '
                f'{named}.  This is the verification pass that runs with `resolve: false` after '
                f'the loop, so these were not cleared by two rounds of rotation followed by a '
                f'minimize.  A doubly-threaded molecule is the usual cause -- one rigid rotation '
                f'about one hinge cannot undo two threadings -- and so is a composition whose '
                f'SAPL is far from where that leaflet sits.  Resolving it here would mean '
                f'starting dynamics on coordinates last changed by an unrelaxed rigid move.')

        # ITERATE.  Clearing one piercing can create another -- a rotated tail has to go
        # somewhere -- and a molecule can be threaded twice, which one rigid rotation about one
        # hinge cannot undo.  So resolve, re-check, and repeat while progress is being made.
        for _cycle in range(_MAX_RESOLVE_CYCLES):
            rotatable = [p for p in npiercings if self._rotatable(p)]
            if not rotatable:
                break
            ess = 's' if len(rotatable) > 1 else ''
            logger.info(f'{len(rotatable)} ring piercing{ess} to resolve by rotation '
                        f'(pass {_cycle + 1}):')
            self._report(rotatable, logger.info)
            fixed_pdb, unresolved = self._resolve_by_rotation(state, rotatable, cutoff,
                                                              segtypes, max_ring_size)
            n_rotations += len(rotatable) - len(unresolved)
            for failed in unresolved:
                logger.info(f'  {self._fmt(failed["piercee"])}: no rotation cleared the piercing '
                            f'by {self._fmt(failed["piercer"])}')
            if fixed_pdb is None:
                break   # nothing moved, so another pass would do exactly the same thing
            # a rotation changes coordinates only, not topology: keep the PSF, new PDB
            self.next_basename('ring_check')
            shutil.copy(fixed_pdb, f'{self.basename}.pdb')
            # the winning poses are the registered artifact now; the per-candidate trial PDBs
            # and VMD gen scripts/logs are scratch
            self._cleanup_declash_scratch()
            self.register(dict(psf=state.psf, pdb=PDBFileArtifact(self.basename),
                               xsc=state.xsc), key='state', artifact_type=StateArtifacts)
            state = self.get_current_artifact('state')
            before = len(npiercings)
            npiercings = run_check(state)
            if not npiercings:
                break
            if len(npiercings) >= before:
                logger.info(f'  {len(npiercings)} piercing(s) remain and that pass made no net '
                            f'progress; stopping the rotation loop')
                break

        # Record the decision rather than leaving it to a log line: a caller that wants to know
        # whether this pass MOVED anything -- and therefore whether a minimize is owed before any
        # dynamics -- must not have to parse INFO output for it.
        self.record_outcome(piercings_initial=n_initial, rotations=n_rotations,
                            piercings_remaining=len(npiercings), resolve=True)

        if not npiercings:
            self.result = 0
            return self.result

        # Whatever remains was not rotation-resolvable.  Lipid rings can be deleted; protein
        # or glycan rings that no rotation cleared are fatal.
        # Anything involving the INSERTED system -- protein or glycan, on EITHER side -- is
        # fatal, and the message says where it came from.  In a membrane build such a piercing
        # can only have been made at the embed, and it cannot survive a working clash test:
        # for a bond to thread a ring it must cross the ring plane inside it, which puts a bond
        # atom 1.60 A from the nearest ring atom for a 6-ring (1.77 for a 7-ring), against the
        # 2.89 A C...C deletion threshold.  Evading that needs a ring radius of 2.79 A, about
        # seventeen atoms; `max_ring_size` is 7.  So the piercing is a report about the clash
        # test, not a thing to resolve -- and resolving it quietly is exactly how a clash test
        # that deleted ZERO of 465 overlapping lipids went unnoticed for three days.
        #
        # A pierced glycan ring additionally has no remedy at all: it cannot be deleted without
        # breaking the glycan tree, which is why the embed clash selection must cover glycans
        # (`noh` over the whole inserted molecule, NOT `protein and noh`).
        def _involves_inserted(p):
            return 'protein' in (p['piercee']['segtype'], p['piercer']['segtype']) \
                or 'glycan' in (p['piercee']['segtype'], p['piercer']['segtype'])

        fatal = [p for p in npiercings if _involves_inserted(p)]
        other_piercings = [p for p in npiercings if not _involves_inserted(p)]

        if fatal:
            self._report(fatal, logger.error)
            named = '; '.join(f'{self._fmt(p["piercee"])} pierced by {self._fmt(p["piercer"])}'
                              for p in fatal)
            raise PestiferBuildError(
                f'{len(fatal)} ring piercing(s) involving protein or glycan: {named}.  A bond '
                f'cannot thread a ring without coming ~1.6 A from a ring atom, far inside the '
                f'embed clash threshold (2.89 A for C...C), so a working protein/glycan clash '
                f'deletion makes this impossible -- check that `bilayer_embed` is actually '
                f'deleting overlapping lipids and solvent (it reports a count), and that its '
                f'heavy-atom selection covers glycans as well as protein.  Rotating this away '
                f'would hide the defect that caused it.')

        if not other_piercings:
            self.result = 0
            return self.result

        ess = 's' if len(other_piercings) > 1 else ''
        if delete_these == "none":
            self._report(other_piercings, logger.error)
            raise PestiferBuildError(
                f'{len(other_piercings)} pierced-ring configuration{ess} found; '
                f'ring_check is configured with delete=none — resolve these manually or change the delete setting'
            )
        # Collect the residues to delete.  A lipid tail threaded through a sterol ring
        # contorts *both* lipids, so a single deletion can leave the other one in a
        # high-energy pose that RATTLE-fails at the first dynamics step; ``both`` removes the
        # piercee and (when it is also a lipid) the piercer, which is the robust choice for a
        # membrane.  The piercer is only deleted if it is a lipid -- never a glycan/protein
        # bond, which cannot be dropped without breaking the molecule.
        to_delete = set()
        for r in other_piercings:
            if delete_these in ('piercee', 'both'):
                to_delete.add((r['piercee']['segname'], str(r['piercee']['resid'])))
            if delete_these in ('piercer', 'both'):
                pr = r['piercer']
                if delete_these == 'piercer' or pr.get('segtype') == 'lipid':
                    to_delete.add((pr['segname'], str(pr['resid'])))
        # deleting a lipid changes the system composition, so report it at INFO (the
        # protein/glycan rotation path is already INFO; a silent lipid delete was misleading)
        logger.info(f'{len(other_piercings)} lipid ring piercing{ess} detected:')
        self._report(other_piercings, logger.info)
        self.next_basename('ring_check')
        pg: PsfgenScripter = self.get_scripter('psfgen')
        pg.newscript(self.basename)
        pg.load_project(state.psf.name, state.pdb.name)
        logger.info(f'Deleting {len(to_delete)} residue(s) (delete={delete_these}) to resolve them:')
        for seg, resid in sorted(to_delete):
            logger.info(f'   deleting segname {seg} residue {resid}')
            pg.addline(f'delatom {seg} {resid}')
        pg.writescript(self.basename)
        pg.runscript()
        self.register(dict(psf=PSFFileArtifact(self.basename), pdb=PDBFileArtifact(self.basename), xsc=state.xsc), key='state', artifact_type=StateArtifacts)
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact)
        self.result = 0
        return self.result

    def _cleanup_declash_scratch(self):
        """Remove the transient files written during rotation resolution: the per-candidate
        trial PDBs (``<taskname>-declash-<seg>-<resid>-chi{1,2}_<deg>.pdb`` and
        ``...-glycan.pdb``) and the VMD gen scripts/logs (``...-gen.tcl``/``...-gen.log``).
        The accepted pose has already been copied to the registered ``<basename>.pdb``
        artifact, so all of these are scratch."""
        for f in glob.glob(f'{self.taskname}-declash-*'):
            try:
                os.remove(f)
            except OSError as e:
                logger.debug(f'could not remove ring_check scratch file {f}: {e}')

    def _resolve_by_rotation(self, state, piercings, cutoff, segtypes, max_ring_size):
        """Clear each rotation-resolvable piercing, one at a time.  The topology is untouched
        (only coordinates move).  Candidates are scored in memory -- a candidate is accepted
        only if it un-threads the ring, and among those the one that introduces the fewest new
        heavy-atom clashes is chosen.

        **Keeps partial progress.**  It used to return on the FIRST piercing it could not
        clear, discarding every rotation already applied: on an example-17 patchB that threw
        away six resolved piercings because the seventh resisted, and reported "8 pierced-ring
        configurations" when two were actually left.  Now it continues, so the caller is handed
        the improved structure and an honest list of what remains.

        Returns ``(fixed_pdb_or_None, unresolved)`` -- ``fixed_pdb`` is ``None`` only when
        nothing at all was written, and ``unresolved`` is the piercings no rotation cleared."""
        psf = state.psf.name
        xsc = state.xsc.name if state.xsc else None
        box = cell_from_xsc(xsc)[0] if xsc else None
        working_pdb = state.pdb.name
        # topology (bond list, ring cycles) is constant across trial rotations -> parse once
        checker = RingChecker(psf, cutoff=cutoff, segtypes=segtypes, max_ring_size=max_ring_size)
        unresolved, resolved_any = [], False
        for p in piercings:
            piercee, piercer = p['piercee'], p['piercer']
            seg, resid, resname = piercee['segname'], piercee['resid'], piercee.get('resname', '?')
            target = (seg, resid)
            base = checker.load_coords(working_pdb)
            out_prefix = f'{self.taskname}-declash-{seg}-{resid}'
            cand = None
            # 1. aromatic ring: rotate the side chain itself
            if piercee['segtype'] == 'protein' and resname in _AROMATIC:
                cand = self._try_sidechain(checker, psf, working_pdb, base, box, xsc,
                                           seg, resid, resname, target, out_prefix)
            # 2. rigid ring (proline / glycan) speared by a glycan: rotate the glycan pendant
            if cand is None and piercer['segtype'] == 'glycan' and piercer.get('bond_serials'):
                cand = self._try_glycan_pendant(checker, working_pdb, base, box, piercer,
                                                target, self._fmt(piercee), out_prefix)
            # 3. glycan ring speared by a protein side chain.  Move the glycan first -- it is the
            #    model-built partner, where the protein side chain generally is not -- by rotating
            #    the branch that carries the pierced ring.  Only if no glycan hinge clears it do we
            #    disturb the protein, rotating the piercing side chain instead.
            if (cand is None and piercee['segtype'] == 'glycan'
                    and piercer['segtype'] == 'protein'):
                cand = self._try_pierced_glycan(checker, working_pdb, base, box, piercee,
                                                target, out_prefix)
                if cand is None:
                    cand = self._try_piercer_sidechain(checker, working_pdb, base, box, piercer,
                                                       target, self._fmt(piercee), out_prefix)
            # 4. lipid ring speared by a lipid bond (a sterol threaded by an acyl tail).
            #    Two moves, cheapest-first by how much they disturb:
            #      a) bend the piercing TAIL about an upstream chain dihedral -- the same
            #         graph-branch rotation the glycan path uses, which is chemistry-agnostic
            #         (its own docstring says so); the moving branch is a few tail carbons.
            #      b) spin the pierced STEROL about its own long axis, in place.
            #    (b) is the fallback because it turns a whole molecule rather than a tail end.
            if (cand is None and piercee['segtype'] == 'lipid'
                    and piercer['segtype'] == 'lipid'):
                if piercer.get('bond_serials'):
                    cand = try_glycan_pendant(checker, working_pdb, base, box, piercer, target,
                                              f'{out_prefix}-lipidtail.pdb', self._fmt(piercee),
                                              mover_kind='lipid tail')
                if cand is None:
                    cand = try_sterol_spin(checker, working_pdb, base, box, piercee, target,
                                           f'{out_prefix}-spin.pdb', self._fmt(piercer))
            # 5. lipid ring speared by a PROTEIN bond -- a sterol threaded by a side chain of an
            #    embedded protein.  Spin the sterol first and disturb the protein only if that
            #    fails, which is case 3's rule applied here: move the model-built partner before
            #    the deposited one, and in a membrane build the lipid is the grid-placed one.
            #    `try_sterol_spin` is honest about its own limit -- it clears a bond passing near
            #    the ring rim and returns None for one threaded through the centre -- so the
            #    fall-through costs a single attempt and the side-chain rotation is the move that
            #    actually pulls the piercing bond out of the ring.
            if (cand is None and piercee['segtype'] == 'lipid'
                    and piercer['segtype'] == 'protein'):
                cand = try_sterol_spin(checker, working_pdb, base, box, piercee, target,
                                       f'{out_prefix}-spin.pdb', self._fmt(piercer))
                if cand is None:
                    cand = self._try_piercer_sidechain(checker, working_pdb, base, box, piercer,
                                                       target, self._fmt(piercee), out_prefix)
            if cand is None:
                unresolved.append(p)
                continue
            working_pdb = cand
            resolved_any = True
        return (working_pdb if resolved_any else None), unresolved

    def _try_sidechain(self, checker, psf, working_pdb, base, box, xsc, seg, resid, resname, target, out_prefix):
        """Rotate the residue's side chain (chi2 then chi1); among the rotations that clear
        the ring, keep the one with the fewest new heavy-atom clashes.  Returns the winning
        PDB or ``None``."""
        best = None  # (clash, pdb, chi, deg)
        for chi, degrees in ((2, _CHI2_DEGREES), (1, _CHI1_DEGREES)):
            prefix = f'{out_prefix}-chi{chi}'
            self._write_sidechain_candidates(psf, working_pdb, seg, resid, chi, degrees, prefix)
            for deg in degrees:
                cand = f'{prefix}_{deg}.pdb'
                if not os.path.exists(cand):
                    continue
                coords = checker.load_coords(cand)
                if checker.check_coords(coords, box, [target]):
                    continue  # still pierced
                moved = np.where(np.abs(coords - base).max(axis=1) > 1e-3)[0]
                # same eligibility rule as the two moves in `ring_resolve`: reject a rotamer
                # that creates a steric overlap rather than ranking it as the least bad
                if checker.clash_count(coords, moved, cutoff=HARD_OVERLAP, box=box):
                    continue
                clash = checker.clash_count(coords, moved, box=box)
                if best is None or clash < best[0]:
                    best = (clash, cand, chi, deg)
                if clash == 0:
                    break
            if best is not None and best[0] == 0:
                break
        if best is None:
            return None
        clash, cand, chi, deg = best
        logger.info(f'  {resname} {seg}-{resid}: chi{chi} rotation of {deg} deg clears the '
                    f'piercing ({clash} new heavy-atom clash(es))')
        return cand

    def _try_glycan_pendant(self, checker, working_pdb, base, box, piercer, target, piercee_label, out_prefix):
        """Rotate the glycan sub-branch carrying the piercing bond about an upstream rotatable
        bond so the piercing bond is pulled out of the ring.  Delegates to the shared
        :func:`~pestifer.psfutil.ring_resolve.try_glycan_pendant` (also used at glycan graft time).
        Returns the winning PDB or ``None``."""
        return try_glycan_pendant(checker, working_pdb, base, box, piercer, target,
                                  f'{out_prefix}-glycan.pdb', piercee_label)

    @staticmethod
    def _sidechain_serials(checker, segname, resid):
        """Serials of one residue's side chain (everything that is not backbone)."""
        return [a.serial for a in checker._atoms
                if a.segname == segname and str(a.resid.resid) == str(resid)
                and a.atomname not in _BACKBONE]

    def _try_pierced_glycan(self, checker, working_pdb, base, box, piercee, target, out_prefix):
        """Swing the *pierced* glycan ring off the bond spearing it, by rotating the glycan
        branch that carries the ring about an upstream glycosidic bond.

        This is the mirror of :meth:`_try_glycan_pendant`, which rotates the branch carrying the
        piercing *bond*.  The same engine does both: it moves whichever graph branch the serials
        it is handed belong to.  A ring is 2-connected, so no bridge can separate its atoms and
        any two of them identify that branch exactly as a bonded pair does."""
        ring_serials = piercee.get('ring_serials') or []
        if len(ring_serials) < 2:
            return None
        as_mover = dict(piercee, bond_serials=list(ring_serials[:2]))
        return try_glycan_pendant(checker, working_pdb, base, box, as_mover, target,
                                  f'{out_prefix}-pierced-glycan.pdb', self._fmt(piercee),
                                  mover_kind='pierced glycan')

    def _try_piercer_sidechain(self, checker, working_pdb, base, box, piercer, target,
                               piercee_label, out_prefix):
        """Rotate the protein side chain that carries the piercing bond, pulling it out of the
        ring.  Restricting the eligible hinges to bonds whose moving branch stays inside the side
        chain leaves exactly that residue's chi bonds, so this works for any residue without a
        per-residue chi table, and the engine's smallest-branch-first ordering tries the
        least-disturbing rotation first."""
        if not piercer.get('bond_serials'):
            return None
        movers = self._sidechain_serials(checker, piercer['segname'], piercer['resid'])
        if not movers:
            return None
        return try_glycan_pendant(checker, working_pdb, base, box, piercer, target,
                                  f'{out_prefix}-sidechain.pdb', piercee_label,
                                  mover_serials=movers, mover_kind='side chain')

    def _write_sidechain_candidates(self, psf, pdb, segname, resid, chi, degrees, out_prefix):
        """Emit one PDB per candidate chi rotation of the residue's side chain
        (``<out_prefix>_<deg>.pdb``), each measured from the input pose -- pure numpy, no VMD."""
        from ..molecule.coordmanip import CoordManipulator
        cm = CoordManipulator(psf, pdb)
        rn = cm.residue_of_segname(segname, resid)
        orig = cm.coords.copy()
        for deg in degrees:
            cm.coords = orig.copy()
            cm.apply_scrot(chi, rn, deg)
            cm.write_pdb(f'{out_prefix}_{deg}.pdb')
