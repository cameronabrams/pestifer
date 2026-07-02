# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`RingCheckTask` class for checking for pierced rings in a molecular structure.
This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and is used to check for pierced rings in a molecular structure.
It identifies configurations where a ring is pierced by a bond, and resolves each according to the pierced ring's segtype:

- a pierced **lipid** ring deletes the offending segment (the piercer or the piercee, chosen by the ``delete`` option; ``none`` reports and stops);
- a pierced **glycan** ring is fatal -- a glycan cannot be deleted without breaking the glycan tree, so the build stops and must be rebuilt (e.g. with a different random seed or more minimization);
- a pierced **protein** side-chain ring (His/Phe/Tyr/Trp/Pro) is resolved by rotating the side chain (chi2, then chi1) until it clears; if no rotation clears it, the build stops with the residue named.

Which segtypes are checked is set by ``segtypes`` (there is no default, so it must be given explicitly).
Usage is described in the :ref:`config_ref tasks ring_check` documentation.
"""
import logging
import os
import shutil

from .basetask import BaseTask
from ..scripters import PsfgenScripter
from ..core.artifacts import *
from ..core.errors import PestiferBuildError
from ..psfutil.psfring import ring_check, RingChecker

logger=logging.getLogger(__name__)

# side-chain dihedral rotations tried, in order, to swing a pierced ring off the
# piercing bond; chi2 (about CB-CG) moves an aromatic ring most directly, chi1 is the
# fallback
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

        # Protein side-chain rings cannot be deleted (they are the structure), so rotate
        # the offending side chain out of the way.  Do this first, then re-check, so any
        # remaining lipid/glycan piercings are handled against the rotated coordinates.
        protein_piercings = [p for p in npiercings if p['piercee']['segtype'] == 'protein']
        if protein_piercings:
            ess = 's' if len(protein_piercings) > 1 else ''
            logger.info(f'{len(protein_piercings)} protein side-chain ring piercing{ess} detected; '
                        f'attempting to rotate the side chain(s) clear:')
            self._report(protein_piercings, logger.info)
            fixed_pdb = self._resolve_protein_piercings(state, protein_piercings, cutoff,
                                                        segtypes, max_ring_size)
            if fixed_pdb is None:
                self._report(protein_piercings, logger.error)
                raise PestiferBuildError(
                    f'{len(protein_piercings)} protein side-chain ring piercing{ess} could not be '
                    f'rotated clear.  Suggested fix: manually rotate the named side-chain dihedral(s) '
                    f'(chi2, then chi1) to swing the ring off the piercing bond, or rebuild with more '
                    f'minimization before ring_check.')
            # a rotation changes coordinates only, not topology: keep the PSF, new PDB
            self.next_basename('ring_check')
            shutil.copy(fixed_pdb, f'{self.basename}.pdb')
            self.register(dict(psf=state.psf, pdb=PDBFileArtifact(self.basename), xsc=state.xsc),
                          key='state', artifact_type=StateArtifacts)
            state = self.get_current_artifact('state')
            npiercings = run_check(state)

        if not npiercings:
            self.result = 0
            return self.result

        glycan_piercings = [p for p in npiercings if p['piercee']['segtype'] == 'glycan']
        remaining_protein = [p for p in npiercings if p['piercee']['segtype'] == 'protein']
        other_piercings   = [p for p in npiercings if p['piercee']['segtype'] not in ('glycan', 'protein')]

        if remaining_protein:
            self._report(remaining_protein, logger.error)
            raise PestiferBuildError(
                f'{len(remaining_protein)} protein side-chain ring piercing(s) remain after rotation; '
                f'resolve manually (rotate the side-chain dihedral) or rebuild')

        if glycan_piercings:
            gess = 's' if len(glycan_piercings) > 1 else ''
            logger.error(f'{len(glycan_piercings)} glycan ring piercing{gess} detected — glycan residues cannot be deleted to resolve:')
            self._report(glycan_piercings, logger.error)
            raise PestiferBuildError(
                f'{len(glycan_piercings)} glycan ring piercing{gess} detected; '
                f'rebuild the system to resolve (e.g. use a different random seed or add more minimization before ring_check)'
            )

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
        self.next_basename('ring_check')
        pg: PsfgenScripter = self.get_scripter('psfgen')
        pg.newscript(self.basename)
        pg.load_project(state.psf.name, state.pdb.name)
        logger.debug(f'Deleting all {delete_these}s from {len(other_piercings)} pierced-ring configuration{ess}')
        for r in other_piercings:
            logger.debug(f'   Deleting segname {r[delete_these]["segname"]} residue {r[delete_these]["resid"]}')
            pg.addline(f'delatom {r[delete_these]["segname"]} {r[delete_these]["resid"]}')
        pg.writescript(self.basename)
        pg.runscript()
        self.register(dict(psf=PSFFileArtifact(self.basename), pdb=PDBFileArtifact(self.basename), xsc=state.xsc), key='state', artifact_type=StateArtifacts)
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact)
        self.result = 0
        return self.result

    def _resolve_protein_piercings(self, state, protein_piercings, cutoff, segtypes, max_ring_size):
        """Rotate each pierced protein side-chain ring out of the way, one residue at a
        time, re-checking after each candidate rotation.  The topology is untouched (only
        coordinates move).  Returns the path to a fixed PDB if every protein piercing
        clears, otherwise ``None``."""
        psf = state.psf.name
        xsc = state.xsc.name if state.xsc else None
        working_pdb = state.pdb.name
        # the PSF (topology, bond list, ring cycles) does not change as we try rotamers,
        # so parse it once and re-check only coordinates for each candidate
        checker = RingChecker(psf, cutoff=cutoff, segtypes=segtypes, max_ring_size=max_ring_size)
        # unique protein residues whose rings are pierced
        residues, seen = [], set()
        for p in protein_piercings:
            side = p['piercee']
            key = (side['segname'], side['resid'])
            if key not in seen:
                seen.add(key)
                residues.append(side)
        for res in residues:
            seg, resid, resname = res['segname'], res['resid'], res.get('resname', '?')
            cleared = False
            for chi, degrees in ((2, _CHI2_DEGREES), (1, _CHI1_DEGREES)):
                prefix = f'{self.taskname}-declash-{seg}-{resid}-chi{chi}'
                self._write_rotamer_candidates(psf, working_pdb, seg, resid, chi, degrees, prefix)
                for deg in degrees:
                    cand = f'{prefix}_{deg}.pdb'
                    if not os.path.exists(cand):
                        continue
                    remaining = checker.check(cand, xsc=xsc, only_piercees=[(seg, resid)])
                    if not remaining:
                        logger.info(f'  {resname} {seg}-{resid}: chi{chi} rotation of {deg} deg clears the piercing')
                        working_pdb = cand
                        cleared = True
                        break
                if cleared:
                    break
            if not cleared:
                logger.error(f'  {resname} {seg}-{resid}: no chi1/chi2 rotation cleared the ring piercing')
                return None
        return working_pdb

    def _write_rotamer_candidates(self, psf, pdb, segname, resid, chi, degrees, out_prefix):
        """Emit one PDB per candidate chi rotation of the given residue's side chain
        (``<out_prefix>_<deg>.pdb``), each measured from the input pose."""
        vm = self.get_scripter('vmd')
        vm.newscript(f'{out_prefix}-gen', packages=['PestiferCRot'])
        vm.addline(f'mol new {psf} waitfor all')
        vm.addline(f'mol addfile {pdb} waitfor all')
        vm.addline('set molid [molinfo top get id]')
        vm.addline(f'set casel [atomselect $molid "segname {segname} and resid {resid} and name CA"]')
        vm.addline('set rn [lindex [$casel get residue] 0]')
        vm.addline('$casel delete')
        vm.addline('set allsel [atomselect $molid all]')
        vm.addline('set orig [$allsel get {x y z}]')
        vm.addline(f'foreach deg {{{" ".join(str(d) for d in degrees)}}} {{')
        vm.addline('  $allsel set {x y z} $orig')
        vm.addline(f'  PestiferCRot::SCrot_chi{chi} $rn $molid $deg')
        vm.addline('  set w [atomselect $molid all]')
        vm.addline(f'  $w writepdb {out_prefix}_$deg.pdb')
        vm.addline('  $w delete')
        vm.addline('}')
        vm.addline('$allsel delete')
        vm.writescript()
        vm.runscript()
