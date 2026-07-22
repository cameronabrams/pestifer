# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`RingCheckTask` class for checking for pierced rings in a molecular structure.
This class is a descendant of the :class:`BaseTask <pestifer.tasks.basetask.BaseTask>` class and is used to check for pierced rings in a molecular structure.
It identifies configurations where a ring is pierced by a bond and resolves each by whichever motion is available:

- an **aromatic** protein side-chain ring (His/Phe/Tyr/Trp) pivots on its own dihedral, so it is swung off the piercing bond by rotating the side chain (chi2, then chi1);
- a **rigid** ring that cannot itself move -- a **proline** side-chain ring (fused to the backbone) or a **glycan** ring -- when speared by a **glycan** bond is cleared by rotating the *glycan pendant* so the piercing bond is pulled out of the ring;
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
from ..psfutil.ring_resolve import try_glycan_pendant
from ..util.util import cell_from_xsc

logger=logging.getLogger(__name__)

# aromatic protein residues whose side-chain ring pivots freely on chi2/chi1 (proline is
# NOT here: its ring is fused into the backbone, so a chi rotation cannot swing it clear)
_AROMATIC = {'HIS', 'HSD', 'HSE', 'HSP', 'PHE', 'TYR', 'TRP'}
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
        return bool(aromatic or pendant)

    def do(self):
        state: StateArtifacts = self.get_current_artifact('state')
        cutoff: float = self.specs.get('cutoff', 3.5)
        segtypes: list = self.specs.get('segtypes', []) or []
        delete_these: str = self.specs.get('delete', 'piercee')
        max_ring_size: int = self.specs.get('max_ring_size', 7)
        if not segtypes:
            logger.warning('ring_check: no segtypes specified (segtypes has no default) — '
                           'nothing checked. Set e.g. segtypes: [lipid, glycan, protein].')
            self.result = 0
            return self.result

        def run_check(st):
            xsc = st.xsc.name if st.xsc else None
            if xsc is None:
                logger.debug('No XSC in pipeline state — ring_check runs in non-periodic (vacuum) mode')
            return ring_check(st.psf.name, st.pdb.name, xsc, cutoff=cutoff,
                              segtypes=segtypes, max_ring_size=max_ring_size)

        npiercings = run_check(state)
        if not npiercings:
            self.result = 0
            return self.result

        # First resolve everything that a rotation can fix (rings that can't be deleted):
        # aromatic side-chain rings by rotating the ring, rigid rings (proline / glycan)
        # speared by a glycan by rotating the glycan pendant.
        rotatable = [p for p in npiercings if self._rotatable(p)]
        if rotatable:
            ess = 's' if len(rotatable) > 1 else ''
            logger.info(f'{len(rotatable)} ring piercing{ess} to resolve by rotation:')
            self._report(rotatable, logger.info)
            fixed_pdb, failed = self._resolve_by_rotation(state, rotatable, cutoff, segtypes, max_ring_size)
            if fixed_pdb is None:
                self._report([failed], logger.error)
                raise PestiferBuildError(
                    f'ring piercing of {self._fmt(failed["piercee"])} by a bond in '
                    f'{self._fmt(failed["piercer"])} could not be cleared by rotation.  Suggested '
                    f'fix: manually rotate the side-chain dihedral (aromatic ring) or the glycan '
                    f'linkage (rigid ring) to un-thread it, or rebuild with more minimization.')
            # a rotation changes coordinates only, not topology: keep the PSF, new PDB
            self.next_basename('ring_check')
            shutil.copy(fixed_pdb, f'{self.basename}.pdb')
            # the winning pose is now the registered artifact; the per-candidate trial PDBs
            # and VMD gen scripts/logs are scratch -- remove them so they don't litter the
            # run directory (they are kept on the failure path above, for debugging)
            self._cleanup_declash_scratch()
            self.register(dict(psf=state.psf, pdb=PDBFileArtifact(self.basename), xsc=state.xsc),
                          key='state', artifact_type=StateArtifacts)
            state = self.get_current_artifact('state')
            npiercings = run_check(state)

        if not npiercings:
            self.result = 0
            return self.result

        # Whatever remains was not rotation-resolvable.  Lipid rings can be deleted; protein
        # or glycan rings that no rotation cleared are fatal.
        glycan_piercings = [p for p in npiercings if p['piercee']['segtype'] == 'glycan']
        protein_piercings = [p for p in npiercings if p['piercee']['segtype'] == 'protein']
        other_piercings = [p for p in npiercings if p['piercee']['segtype'] not in ('glycan', 'protein')]

        if protein_piercings or glycan_piercings:
            fatal = protein_piercings + glycan_piercings
            self._report(fatal, logger.error)
            raise PestiferBuildError(
                f'{len(fatal)} protein/glycan ring piercing(s) could not be resolved by rotation; '
                f'rebuild the system to resolve (e.g. a different random seed or more minimization '
                f'before ring_check)')

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
        heavy-atom clashes is chosen.  Returns ``(fixed_pdb, None)`` if every piercing clears,
        or ``(None, failed_piercing)`` for the first that does not."""
        psf = state.psf.name
        xsc = state.xsc.name if state.xsc else None
        box = cell_from_xsc(xsc)[0] if xsc else None
        working_pdb = state.pdb.name
        # topology (bond list, ring cycles) is constant across trial rotations -> parse once
        checker = RingChecker(psf, cutoff=cutoff, segtypes=segtypes, max_ring_size=max_ring_size)
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
            if cand is None:
                return None, p
            working_pdb = cand
        return working_pdb, None

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
                clash = checker.clash_count(coords, moved)
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
