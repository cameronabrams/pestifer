# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`PsfgenTask` class for handling invocations of psfgen which create a molecule from a base PDB/mmCIF file.

Usage is described in the :ref:`subs_buildtasks_psfgen` documentation.

"""
import logging
import networkx as nx
import numpy as np
import os
import shutil

from copy import deepcopy
from pathlib import Path
from typing import ClassVar


from .basetask import VMDTask
from ..molecule.chainidmanager import ChainIDManager
from ..molecule.molecule import Molecule
from ..core.objmanager import ObjManager
from ..core.artifacts import *
from ..core.errors import PestiferBuildError
from pidibble.pdbparse import PDBParser

from ..objs.cfusion import CfusionList
from ..objs.graft import GraftList
from ..psfutil.psfatom import PSFAtomList
from ..psfutil.psfcontents import PSFContents
from ..scripters import VMDScripter, PsfgenScripter
from ..util.util import write_residue_map

logger = logging.getLogger(__name__)

class PsfgenTask(VMDTask):
    """ 
    A class for handling invocations of psfgen which create a molecule from a base PDB/mmCIF file
    or from a PSF file generated previously by psfgen
    This class is a descendant of the :class:`BaseTask <pestifer.tasks.basetask.BaseTask>` class, and as such it has two parameters.

    Parameters
    ----------
    config_specs : dict
        Configuration specifications for the task.
    """

    _yaml_header: ClassVar[str] = 'psfgen'
    """
    YAML header for the PsfgenTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, SOURCE, STATE, MOLECULE
        # builds a system from fetched source coordinates; running it on top of an existing
        # state discards that state (it rebuilds from scratch)
        return TaskContract(requires=(SOURCE,), provides=(STATE, MOLECULE), discards_state=True)

    def provision(self, packet: dict = {}):
        """
        Provision the task with the standard packet, and initialize molecule dict and base molecule placeholders.
        """
        super().provision(packet=packet)
        self.molecules = {}
        self.base_molecule: Molecule = None

    def do(self):
        """
        Execute the psfgen task.
        This method initializes the task, ingests the base molecule, and runs the psfgen process.
        It also handles any necessary coormods and declashing of loops and glycans based on the task specifications.
        The results of the psfgen process are saved as a PSF/PDB fileset, and the state is updated accordingly.
        """
        logger.debug('ingesting molecule(s)')
        self.ingest_molecules()
        logger.debug(f'base mol num images {self.base_molecule.num_images()}')
        self._check_ion_aliases()
        logger.debug('Running first psfgen')
        self.result = self.psfgen()
        if self.result != 0:
            return self.result
        # we now have a full coordinate set, so we can do coormods
        self.coormods()
        self.declash()
        return self.result

    def coormods(self):
        """
        Perform coordinate modifications based on the specifications provided in the task.
        """
        coormods = self.objmanager.get('coord',{})
        logger.debug(f'psfgen task has {len(coormods)} coormods:')
        logger.debug(';'.join([str(_) for _ in coormods]))

        if coormods:
            logger.debug(f'performing coormods')
            for objtype, objlist in coormods.items():
                if len(objlist) == 0:
                    continue
                self.next_basename(objtype)
                state: StateArtifacts = self.get_current_artifact('state')
                if objtype in ('irotations', 'crotations'):
                    # internal-coordinate (torsion) rotations -- applied in pure numpy, once per
                    # biological-assembly image (chainIDmap remaps the crotation onto each copy).
                    if self._apply_python_crotations(objlist, state) != 0:
                        self.result = 1
                        return
                    self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf,
                                       xsc=state.xsc), key='state', artifact_type=StateArtifacts)
                    continue
                if objtype == 'orient':
                    # principal-axis reorientation of the whole molecule -- pure numpy
                    if self._apply_python_orient(objlist, state) != 0:
                        self.result = 1
                        return
                    self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf,
                                       xsc=state.xsc), key='state', artifact_type=StateArtifacts)
                    continue
                vm: VMDScripter = self.scripters['vmd']
                vm.newscript(self.basename, packages=[])
                vm.load_psf_pdb(state.psf.name, state.pdb.name, new_molid_varname='mCM')
                match objtype:
                    case 'rottrans':
                        vm.write_rottranslist(objlist)
                vm.write_pdb(self.basename, 'mCM')
                vm.writescript()
                vm.runscript()
                self.register(dict(
                    pdb=PDBFileArtifact(self.basename),
                    psf=state.psf,
                    xsc=state.xsc), key='state', artifact_type=StateArtifacts)
                self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
                self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)

    def _apply_python_crotations(self, objlist, state: StateArtifacts) -> int:
        """Apply irotations/crotations to the built psf/pdb in numpy, once per assembly image, and
        write ``{basename}.pdb``.  Returns 0 on success, 1 if a rotation could not be applied."""
        from ..molecule.coordmanip import CoordManipulator, CoordManipulateError
        cm = CoordManipulator(state.psf.name, state.pdb.name)
        try:
            for transform in self.base_molecule.active_biological_assembly.transforms.data:
                for crot in objlist:
                    cm.apply_crot(crot, chainIDmap=transform.chainIDmap)
        except CoordManipulateError as e:
            logger.error(f'psfgen crotations: {e}')
            return 1
        cm.write_pdb(f'{self.basename}.pdb')
        return 0

    def _apply_python_orient(self, objlist, state: StateArtifacts) -> int:
        """Apply principal-axis orient(s) to the built psf/pdb in numpy, writing ``{basename}.pdb``."""
        from ..molecule.coordmanip import CoordManipulator, CoordManipulateError
        cm = CoordManipulator(state.psf.name, state.pdb.name)
        try:
            for orient in objlist:
                cm.apply_orient(orient)
        except CoordManipulateError as e:
            logger.error(f'psfgen orient: {e}')
            return 1
        cm.write_pdb(f'{self.basename}.pdb')
        return 0

    def declash(self):
        """
        Manages the declashing of protein loops and glycans
        """
        self.min_loop_length = self.specs['source'].get('sequence',{}).get('loops',{}).get('min_loop_length',0)
        self.declash_counts = self.base_molecule.loop_counts(min_loop_length=self.min_loop_length)
        num_images = self.base_molecule.num_images()
        for segtype in ['protein','nucleicacid']:
            self.declash_counts[segtype] *= num_images
        self.declash_counts['glycan'] = self.base_molecule.nglycans() * num_images
        for segtype,speckey in zip(['protein','nucleicacid','glycan'],['loops','loops','glycans']):
            if self.declash_counts[segtype]>0 and self.specs['source']['sequence'][speckey]['declash']['maxcycles']>0:
                logger.debug(f'Declashing {self.declash_counts[segtype]} {segtype} segments')
                self.declash_segtype(self.specs['source']['sequence'][speckey]['declash'], segtype=segtype)
        # Model built terminal tails (opted in via build_zero_occupancy_*/include_terminal_loops)
        # with the free-tail modeler, so a tail gets the same realistic, clash-filtered backbone the
        # interior CCD closer gives a loop instead of guesscoord's extended arm.
        self.model_terminal_tails()
        # A grafted glycan bond can thread through a protein/glycan ring with no atomic clash, so
        # the clash-count declash above cannot see it.  Detect such topological piercings and
        # rotate the glycan pendant out of the ring here, before the structure is written --
        # ideally making a downstream ring_check task unnecessary.
        self.resolve_glycan_piercings()

    def model_terminal_tails(self):
        """
        Model built protein terminal tails with the free-tail modeler -- the terminal analogue
        of the interior-loop CCD closer.

        A tail opted in via ``build_zero_occupancy_[NC]_termini`` / ``include_terminal_loops`` is
        emitted by psfgen as an arbitrary extended ``guesscoord`` arm. Here each such tail's
        backbone is replaced by a Ramachandran-seeded, clash-filtered conformation (rigidly
        anchored at its resolved junction; no closure, since a tail has a free end), so terminal
        residues are modeled to the same standard as interior loops. Best-effort and
        coordinate-only (the PSF is unchanged); a downstream ``minimize`` relaxes residual soft
        overlaps, exactly as for interior loops. Toggle ``loops.declash.model_tails`` (default on).
        """
        specs = self.specs['source']['sequence']['loops']['declash']
        if not specs.get('model_tails', True):
            return
        tails = self.base_molecule.protein_terminal_tails()
        if not tails:
            return
        import numpy as np
        from ..psfutil.tail_model import model_one_tail
        from ..psfutil.loop_ccd import loop_atoms_from_pdb
        from ..util.coord import pdb_replace_coords
        self.next_basename('model-tails')
        state: StateArtifacts = self.get_current_artifact('state')
        src_pdb = state.pdb.name
        ensemble = int(specs.get('tail_ensemble', 10))
        # Every model-built residue (interior loops + all tails) holds throwaway guesscoord
        # coords at this point, so exclude all of them from every tail's clash environment --
        # only the resolved structure is a real steric backdrop.
        modeled_specs = [(g['segname'], g['loop_resids']) for g in self.base_molecule.protein_loop_gaps()]
        modeled_specs += [(t['segname'], t['tail_resids']) for t in tails]
        all_modeled_serials = set()
        for segname, resids in modeled_specs:
            _o, _c, ser = loop_atoms_from_pdb(src_pdb, resids, segname=segname)
            all_modeled_serials.update(ser)
        all_modeled_serials = sorted(all_modeled_serials)
        new_coords = {}
        placed_tails = []      # heavy coords of tails already modeled this pass
        rotation_report = []   # (tag, segname, rotation rows) per tail, for the provenance file
        tail_serials = set()   # all model-built tail atom serials (movers for piercing resolution)
        for k, t in enumerate(tails):
            extra_env = np.vstack(placed_tails) if placed_tails else None
            res = model_one_tail(src_pdb, t['segname'], t['tail_resids'], t['anchor_resid'],
                                 t['end'], all_modeled_serials,
                                 seed=self._DECLASH_SEED + k, ensemble=ensemble, extra_env=extra_env)
            placed_tails.append(res['heavy'])
            rep = res['rep']
            tag = (f"{t['segname']}:{t['tail_resids'][0]}-{t['tail_resids'][-1]} "
                   f"({t['end']}-term, len {len(t['tail_resids'])})")
            rotation_report.append((tag, t['segname'], res.get('rotations', [])))
            deep = rep['n_deep'] + rep['n_env_deep']
            soft = rep['n_soft'] + rep['n_env_soft']
            if deep:
                logger.info(f'psfgen/model-tails: tail {tag} modeled with {deep} deep + {soft} soft '
                            f'overlaps (worst {rep["worst"]:.2f} A) -- a downstream minimize relaxes these.')
            else:
                logger.debug(f'psfgen/model-tails: tail {tag} modeled clash-free ({soft} soft).')
            for j, serial in enumerate(res['serials']):
                new_coords[serial] = res['modeled'][j]
            tail_serials.update(res['serials'])
        out_pdb = f'{self.basename}.pdb'
        serial_list = sorted(new_coords)
        coords_arr = np.array([new_coords[s] for s in serial_list])
        row_of_serial = {s: i for i, s in enumerate(serial_list)}
        pdb_replace_coords(src_pdb, out_pdb, coords_arr, row_of_serial)
        self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf,
                           xsc=getattr(state, 'xsc', None)), key='state', artifact_type=StateArtifacts)
        self.write_rotation_report(rotation_report, kind='terminal tail modeling')
        logger.info(f'psfgen: modeled {len(tails)} terminal tail(s) -> {self.basename}.pdb')
        if specs.get('check_piercings', True):
            self.resolve_tail_piercings(sorted(tail_serials))

    def resolve_tail_piercings(self, tail_serials):
        """Detect and clear any ring piercing caused by a model-built terminal-tail backbone bond,
        by rotating the offending tail sub-branch out of the ring.

        A threading has no atomic clash, so the tail modeler's declash cannot see it; this scan
        catches the rare case where a modeled tail spears a proline/aromatic/sugar ring. A tail has
        a free end, so swinging its distal portion about a backbone bridge is unconstrained. Shares
        the glycan-graft resolution engine (:func:`~pestifer.psfutil.ring_resolve.resolve_model_built_piercings`);
        best-effort and non-fatal (toggle ``loops.declash.check_piercings``).
        """
        if not tail_serials:
            return
        from ..psfutil.ring_resolve import resolve_model_built_piercings
        state: StateArtifacts = self.get_current_artifact('state')
        xsc = state.xsc.name if getattr(state, 'xsc', None) else None
        self.next_basename('tail-piercing-check')
        fixed_pdb, unresolved = resolve_model_built_piercings(
            state.psf.name, state.pdb.name, xsc, self.basename, tail_serials,
            segtypes=['protein', 'glycan'], mover_kind='tail')
        if fixed_pdb is not None:
            shutil.copy(fixed_pdb, f'{self.basename}.pdb')
            self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf,
                               xsc=getattr(state, 'xsc', None)),
                          key='state', artifact_type=StateArtifacts)
            logger.info(f'Resolved terminal-tail ring-piercing(s) -> {self.basename}.pdb')
        if unresolved:
            logger.warning(f'{len(unresolved)} terminal-tail ring-piercing(s) could not be cleared '
                           f'by rotation; add a ring_check task or rebuild:')
            for p in unresolved:
                ee, er = p['piercee'], p['piercer']
                logger.warning(f'  ring of {ee.get("resname", "?")} {ee["segname"]}-{ee["resid"]} '
                               f'pierced by a bond in tail {er["segname"]}-{er["resid"]}')

    def resolve_glycan_piercings(self):
        """Detect ring piercings caused by grafted/model-built glycan bonds and clear each by
        rotating the offending glycan sub-branch out of the ring (see
        :func:`~pestifer.psfutil.ring_resolve.resolve_glycan_piercings`).

        Runs whenever the base molecule carries glycans, independent of the clash-count declash
        cycles (a piercing has no atomic clash).  Disabled with
        ``glycans.declash.check_piercings: false``.  Best-effort and non-fatal: a piercing no
        rotation clears is reported as a warning, leaving a downstream ``ring_check`` (if any) to
        attempt further resolution.
        """
        declash_specs = self.specs.get('source', {}).get('sequence', {}).get('glycans', {}).get('declash', {})
        if not declash_specs.get('check_piercings', True):
            return
        if self.declash_counts.get('glycan', 0) <= 0:
            return
        from ..psfutil.ring_resolve import resolve_glycan_piercings as _resolve
        state: StateArtifacts = self.get_current_artifact('state')
        xsc = state.xsc.name if getattr(state, 'xsc', None) else None
        self.next_basename('glycan-piercing-check')
        logger.debug('Checking grafted glycans for ring piercings')
        fixed_pdb, unresolved = _resolve(state.psf.name, state.pdb.name, xsc, self.basename,
                                         segtypes=['protein', 'glycan'])
        if fixed_pdb is not None:
            shutil.copy(fixed_pdb, f'{self.basename}.pdb')
            self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf,
                               xsc=getattr(state, 'xsc', None)),
                          key='state', artifact_type=StateArtifacts)
            logger.info(f'Resolved glycan ring-piercing(s) at graft time -> {self.basename}.pdb')
        if unresolved:
            logger.warning(f'{len(unresolved)} ring-piercing(s) could not be cleared by glycan '
                           f'rotation at graft time; add a ring_check task to attempt further '
                           f'resolution (e.g. aromatic side-chain rotation) or rebuild:')
            for p in unresolved:
                ee, er = p['piercee'], p['piercer']
                logger.warning(f'  ring of {ee.get("resname","?")} {ee["segname"]}-{ee["resid"]} '
                               f'pierced by a bond in {er.get("resname","?")} '
                               f'{er["segname"]}-{er["resid"]} ({er["segtype"]})')

    def resi_topologies(self):
        """
        Collect the topology files that are needed for the residues in the base molecule.
        """
        resis = set([x.resname for x in self.base_molecule.asymmetric_unit.residues.data])
        logger.debug(f'Residue names in base molecule: {resis}')
        CC = self.resource_manager.charmmff_content
        new_topfiles = set()
        for resname in resis:
            topfile = CC.get_topfile_of_resname(resname)
            if topfile:
                new_topfiles.add(topfile)
            else:
                raise PestiferBuildError(f'No topology file found for residue name {resname}')
        return list(new_topfiles)

    def patch_topologies(self):
        """
        Collect the topology files that are needed for the patches and links in the base molecule.
        This method retrieves the CHARMMFF topology files associated with the patches and links defined 
        in the base molecule's object manager.
        """
        objmanager = self.base_molecule.objmanager
        seqmods = objmanager.get('seq', {})
        patches = seqmods.get('patches', [])
        CC = self.resource_manager.charmmff_content
        new_topfiles = set()
        # logger.debug(f'New topologies: {new_topfiles}')
        for patch in patches:
            topfile = CC.get_topfile_of_resname(patch.patchname)
            if topfile:
                new_topfiles.add(topfile)
        # logger.debug(f'New topologies: {new_topfiles}')
        topomods = objmanager.get('topol', {})
        links = topomods.get('links', [])
        for link in links:
            topfile = CC.get_topfile_of_resname(link.patchname)
            if topfile:
                new_topfiles.add(topfile)
        # logger.debug(f'New topologies: {new_topfiles}')
        return list(new_topfiles)

    def psfgen(self):
        """
        Run the psfgen process to generate a PSF and PDB fileset from the base molecule.
        """
        self.next_basename('build')
        pg: PsfgenScripter = self.scripters['psfgen']
        required_topology_files = list(set(self.patch_topologies() + self.resi_topologies()))
        pg.newscript(self.basename, packages=['PestiferCRot'], additional_topologies=required_topology_files)
        pg.set_molecule(self.base_molecule, altcoords=self.specs.get('source', {}).get('altcoords', None))
        pg.describe_molecule(self.base_molecule)
        pg.writescript(self.basename)
        result = pg.runscript(keep_tempfiles=True)
        if result != 0:
            return result
        # register PSF, PDB, log, and all charmmff files in the pipeline context
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact) 
        self.register(dict(
            pdb=PDBFileArtifact(self.basename, pytestable=True), 
            psf=PSFFileArtifact(self.basename, pytestable=True)), key='state', artifact_type=StateArtifacts)
        self.register([CharmmffTopFileArtifact(x) for x in pg.topologies if x.endswith('.rtf')], key='charmmff_topfiles', artifact_type=CharmmffTopFileArtifacts)
        self.register([CharmmffStreamFileArtifact(x) for x in pg.topologies if x.endswith('.str')], key='charmmff_streamfiles', artifact_type=CharmmffStreamFileArtifacts)
        self.register([PDBFileArtifact(x) for x in pg.F if x.endswith('.pdb')], key='psfgen_temp_pdbs', artifact_type=PDBFileArtifactList)
        self.pipeline.show_artifacts(header='Artifacts after psfgen')
        self.strip_remarks()
        self._assert_coordinates_set()
        return 0

    # common monatomic ions whose PDB residue name differs from (and, as with calcium, can collide
    # with) the CHARMM residue name.  A bare (monatomic) residue with one of these PDB names is
    # unambiguously the ion, but psfgen would build the same-named CHARMM/CGenFF residue instead --
    # so pestifer stops early and tells the user the alias to add.  Mapping: PDB resname ->
    # (CHARMM resname, element).  For a monatomic ion the atom name equals the residue name on both
    # sides, so the atom alias is "<CHARMM> <PDB> <CHARMM>".
    _ION_PDB_TO_CHARMM = {
        'CA':  ('CAL', 'calcium'),
        'NA':  ('SOD', 'sodium'),
        'K':   ('POT', 'potassium'),
        'CL':  ('CLA', 'chloride'),
        'ZN':  ('ZN2', 'zinc'),
        'CS':  ('CES', 'cesium'),
        'RB':  ('RUB', 'rubidium'),
        'LI':  ('LIT', 'lithium'),
        'CD':  ('CD2', 'cadmium'),
    }

    def _check_ion_aliases(self):
        """Hard-error *before* running psfgen if the input contains a monatomic ion whose PDB
        residue name must be aliased to a differently-named CHARMM residue (e.g. calcium ``CA`` ->
        ``CAL``) and no such residue alias is configured.

        Left un-aliased, psfgen builds the same-named CHARMM/CGenFF residue (a multi-atom organic
        molecule, for ``CA``) and cannot place its extra atoms -- which the origin-coordinate guard
        would later catch, but only after a wasted build.  This stops early with the exact aliases
        to add, and never silently reinterprets the user's input.
        """
        pg: PsfgenScripter = self.scripters['psfgen']
        residue_aliases = (getattr(pg, 'psfgen_config', {}) or {}).get('aliases', {}).get('residue', []) or []
        already_aliased = {entry.split()[0] for entry in residue_aliases if entry.split()}
        problems = {}  # pdb_resname -> [residue labels]
        for res in self.base_molecule.asymmetric_unit.residues.data:
            rn = res.resname
            if rn in self._ION_PDB_TO_CHARMM and rn not in already_aliased and len(res.atoms) == 1:
                problems.setdefault(rn, []).append(f'{res.chainID}{res.resid.resid}')
        if not problems:
            return
        described = '; '.join(
            f'{rn} ({self._ION_PDB_TO_CHARMM[rn][1]}) at {", ".join(labels)}'
            for rn, labels in sorted(problems.items()))
        res_aliases = ', '.join(f'"{rn} {self._ION_PDB_TO_CHARMM[rn][0]}"' for rn in sorted(problems))
        atom_aliases = ', '.join(f'"{self._ION_PDB_TO_CHARMM[rn][0]} {rn} {self._ION_PDB_TO_CHARMM[rn][0]}"'
                                 for rn in sorted(problems))
        raise PestiferBuildError(
            f'The input contains monatomic ion residue(s) whose PDB name must be aliased to the '
            f'CHARMM ion residue before psfgen (otherwise psfgen builds the wrong, same-named '
            f'residue and cannot place its atoms): {described}.\n'
            f'Add to your config:\n'
            f'    psfgen:\n'
            f'      aliases:\n'
            f'        residue: [{res_aliases}]\n'
            f'        atom: [{atom_aliases}]')

    def _assert_coordinates_set(self):
        """Hard-error if the psfgen output PDB contains atoms left at the origin (0,0,0).

        psfgen writes any atom whose coordinates it could neither read from the input nor guess at
        exactly ``0.000 0.000 0.000``.  This most often happens when a PDB residue *name* collides
        with a different, multi-atom CHARMM/CGenFF residue of the same name -- classically a
        monatomic metal ion (e.g. calcium named ``CA``, which matches a CGenFF organic residue),
        so psfgen builds the organic residue's extra atoms and has no coordinates for them.  Left
        unchecked these origin atoms make the following minimization blow up, and pestifer would
        otherwise treat the psfgen output as a successful build.
        """
        state: StateArtifacts = self.get_current_artifact('state')
        pdb = state.pdb.name
        if not pdb or not os.path.exists(pdb):
            return
        residues = {}
        with open(pdb) as f:
            for line in f:
                if line[:6] in ('ATOM  ', 'HETATM') and line[30:54] == '   0.000   0.000   0.000':
                    resname, chain, resid = line[17:21].strip(), line[21].strip(), line[22:26].strip()
                    residues[(resname, chain, resid)] = residues.get((resname, chain, resid), 0) + 1
        if not residues:
            return
        n_atoms = sum(residues.values())
        reslist = ', '.join(f'{rn} {ch}{rid} ({n} atom{"s" if n != 1 else ""})'
                            for (rn, ch, rid), n in sorted(residues.items()))
        raise PestiferBuildError(
            f'psfgen left {n_atoms} atom(s) at the origin (0,0,0) in {pdb}: their coordinates could '
            f'not be read from the input or guessed. Affected residue(s): {reslist}.\n'
            f'This usually means a PDB residue name collides with a different, multi-atom '
            f'CHARMM/CGenFF residue of the same name -- most often a monatomic metal ion (e.g. '
            f'calcium named CA, which matches a CGenFF organic residue). Alias it to the CHARMM ion '
            f'residue in your config; for calcium that is:\n'
            f'    psfgen:\n'
            f'      aliases:\n'
            f'        residue: ["CA CAL"]\n'
            f'        atom: ["CAL CA CAL"]')

    def strip_remarks(self):
        """
        Strip REMARK lines from the PDB file generated by psfgen.
        This method removes any REMARK lines from the PDB file to ensure that it contains only the relevant atomic coordinates and structure information.
        """
        state: StateArtifacts = self.get_current_artifact('state')
        pdb = state.pdb.name
        if not pdb:
            logger.warning('No PDB file found to strip remarks from.')
            return
        if not os.path.exists(pdb):
            logger.warning(f'PDB file {pdb} does not exist, cannot strip remarks.')
            return
        logger.debug(f'Stripping REMARK lines from {pdb}.')
        with open(pdb, 'r') as infile:
            lines = infile.readlines()
        with open(pdb, 'w') as outfile:
            for line in lines:
                if not line.startswith('REMARK'):
                    outfile.write(line)
        self.register(dict(
            pdb=PDBFileArtifact(self.basename, pytestable=True), 
            psf=state.psf), key='state', artifact_type=StateArtifacts)

    def declash_segtype(self, specs: dict, segtype='protein'):
        """
        Dispatch declashing of a segment type (protein gap loops, nucleic-acid loops, or glycans)
        to the corresponding pure-numpy declasher.
        """
        specs_to_declashers = specs#['source']['sequence']['declash']
        if segtype == 'protein':
            self.declash_protein_loops(specs_to_declashers)
        elif segtype == 'nucleicacid':
            self.declash_na_loops(specs_to_declashers)
        elif segtype == 'glycan':
            self.declash_glycans(specs_to_declashers)
        else:
            raise PestiferBuildError(f'Unknown segment type {segtype} for declashing')

    def declash_protein_loops(self, specs: dict):
        """
        Declash model-built protein gap loops in pure numpy (no VMD): greedily wiggle each loop's
        backbone phi angles to pull the loop out of steric overlap with the surrounding structure.
        """
        cycles = specs['maxcycles']
        if not cycles:
            return
        self.next_basename('declash-loops')
        state: StateArtifacts = self.get_current_artifact('state')
        # When the free-tail modeler is on it owns every terminal tail (N- and C-terminal), so the
        # phi-wiggle handles interior loops only and does not also wiggle a C-terminal tail.
        include_c = specs['include_C_termini'] and not specs.get('model_tails', True)
        loops = self._enumerate_protein_loops(include_c_termini=include_c)
        logger.debug(f'Declashing {len(loops)} protein loops (image copies); maxcycles {cycles}')
        self._run_loop_declash(state, loops, cycles)
        self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf),
                      key='state', artifact_type=StateArtifacts)

    def _enumerate_protein_loops(self, include_c_termini=False):
        """Return ``[(segname, [resids]), ...]`` for each model-built MISSING protein loop, once per
        biological-assembly image (segname is the image's chain label)."""
        mol = self.base_molecule
        ba = mol.active_biological_assembly
        au = mol.asymmetric_unit
        min_length = self.min_loop_length
        out = []
        for S in (s for s in au.segments.data if s.segtype == 'protein'):
            n_subsegs = len(S.subsegments)
            for b in S.subsegments.data:
                is_c_terminus = (S.subsegments.index(b) == (n_subsegs - 1))
                is_processible = b.state == 'MISSING' and b.num_items() >= min_length
                if is_processible and (not include_c_termini) and is_c_terminus:
                    is_processible = False
                if not is_processible:
                    continue
                resids = [r.resid.resseqnum for r in S.residues.data[b.bounds[0]:b.bounds[1] + 1]]
                for transform in ba.transforms:
                    out.append((transform.chainIDmap.get(S.segname, S.segname), resids))
        return out

    def _run_loop_declash(self, state, loops, cycles):
        """Load the current psf/pdb and run the numpy phi-wiggle loop declasher over each
        ``(segname, resids)`` loop against its static (non-loop) environment; write ``{basename}.pdb``."""
        from ..molecule.coordmanip import CoordManipulator, CoordManipulateError
        from ..psfutil.declash import declash_loop
        cm = CoordManipulator(state.psf.name, state.pdb.name)
        coords = cm.coords
        heavy = cm._mass > 1.1
        ridx = cm._residue_index()
        names = cm._names
        seg = np.array([a.segname for a in cm.atoms.data])
        rng = np.random.default_rng(self._DECLASH_SEED)
        not_nhn = ~np.isin(names, ['N', 'HN'])
        for segname, resids in loops:
            seg_mask = seg == segname
            loop_mask = seg_mask & np.isin([a.resid.resseqnum for a in cm.atoms.data], resids)
            env_heavy = np.nonzero(heavy & ~loop_mask)[0]
            try:
                r_end = cm.residue_of_segname(segname, resids[-1])
            except CoordManipulateError:
                continue
            residues = []
            for ri in resids:
                r0 = cm.residue_of_segname(segname, ri)
                n_i = np.nonzero(seg_mask & (ridx == r0) & (names == 'N'))[0]
                ca_i = np.nonzero(seg_mask & (ridx == r0) & (names == 'CA'))[0]
                if not n_i.size or not ca_i.size:
                    continue
                movers = np.nonzero(seg_mask & (((ridx == r0) & not_nhn) | ((ridx > r0) & (ridx <= r_end))))[0]
                frag_heavy = np.nonzero(seg_mask & heavy & (ridx >= r0) & (ridx <= r_end))[0]
                residues.append({'pivot': int(n_i[0]), 'axis_to': int(ca_i[0]),
                                 'movers': movers, 'frag_heavy': frag_heavy})
            declash_loop(coords, residues, env_heavy, cycles, 2.0, rng)
        cm.coords = coords
        cm.write_pdb(f'{self.basename}.pdb')

    def declash_na_loops(self, specs):
        """
        Declash nucleic-acid loops in pure numpy (no VMD): rotate each model-built NA loop's
        rotatable bonds (base groups) out of steric overlap with the rest of the structure.
        """
        cycles = specs['maxcycles']
        clashdist = specs['clashdist']
        self.next_basename('declash-na-loops')
        state: StateArtifacts = self.get_current_artifact('state')
        loops = self._enumerate_na_declash(state.psf.name)
        logger.debug(f'Declashing {len(loops)} nucleic acid loops; clashdist {clashdist}; maxcycles {cycles}')
        self._run_pendant_declash(state, loops, cycles, clashdist)
        self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf),
                      key='state', artifact_type=StateArtifacts)

    def _enumerate_na_declash(self, psf, **options):
        """Return ``[(indices, bonds, movers), ...]`` for each model-built nucleic-acid gap loop."""
        au = self.base_molecule.asymmetric_unit
        struct = PSFContents(psf, parse_topology=['bonds'])
        na_atoms: PSFAtomList = struct.atoms.get(lambda x: x.segtype == 'nucleicacid')
        if not na_atoms:
            return []
        min_length = self.min_loop_length
        include_c_termini = options.get('include_c_termini', False)
        out = []
        for S in (s for s in au.segments.data if s.segtype == 'nucleicacid'):
            n_subsegs = len(S.subsegments)
            for b in S.subsegments.data:
                lr_resid = S.residues[b.bounds[0]].resid
                rr_resid = S.residues[b.bounds[1]].resid
                is_c_terminus = (S.subsegments.index(b) == (n_subsegs - 1))
                is_processible = b.state == 'MISSING' and b.num_items() >= min_length
                if is_processible and (not include_c_termini) and is_c_terminus:
                    is_processible = False
                if not is_processible:
                    continue
                loop_atoms = PSFAtomList([x for x in na_atoms.data if x.segname == S.segname
                                          and lr_resid <= x.resid <= rr_resid])
                G = [loop_atoms.graph().subgraph(c).copy()
                     for c in nx.connected_components(loop_atoms.graph())]
                assert len(G) == 1, 'NA loop has more than one connected component'
                g = G[0]
                serials = {x.serial for x in g}
                for at in g:                          # roots = atoms bonded outside the loop (both ends)
                    if any(ls not in serials for ls in (x.serial for x in at.ligands)):
                        at.is_root = True
                out.append((np.array([x.serial - 1 for x in g], dtype=int),) + self._rotatable_bonds(g))
        return out

    def declash_glycans(self, specs):
        """
        Declash glycans in the base molecule in pure numpy (no VMD): rotate each glycan's rotatable
        bonds to pull the sugar chains out of steric overlap with the rest of the structure before
        minimization.  Greedy MC scored by a heavy-atom contact count, seeded for reproducibility.
        """
        mol: Molecule = self.base_molecule
        cycles: int = specs['maxcycles']
        clashdist: float = specs['clashdist']
        if not mol.nglycans() or not cycles:
            logger.debug(f'Glycan declashing is intentionally not done.')
            return
        self.next_basename('declash-glycans')
        state: StateArtifacts = self.get_current_artifact('state')
        glycans = self._enumerate_glycan_declash(state.psf.name)
        logger.debug(f'Declashing {len(glycans)} glycans; clashdist {clashdist}; maxcycles {cycles}')
        self._run_pendant_declash(state, glycans, cycles, clashdist)
        self.register(dict(pdb=PDBFileArtifact(self.basename), psf=state.psf),
                      key='state', artifact_type=StateArtifacts)

    def _enumerate_glycan_declash(self, psf):
        """Return ``[(indices, bonds, movers), ...]`` for each glycan (0-based atom row indices):
        the pendant atoms, its rotatable bridging bonds, and the mover atoms per bond -- the
        declash inputs formerly emitted as Tcl."""
        struct = PSFContents(psf, parse_topology=['bonds'])
        glycanatoms = struct.atoms.get(lambda x: x.segtype == 'glycan')
        glycangraph = glycanatoms.graph()
        G = [glycangraph.subgraph(c).copy() for c in nx.connected_components(glycangraph)]
        out = []
        for g in G:
            serials = {x.serial for x in g}
            rooted = False
            for at in g:                              # root = the atom bonded outside the glycan
                for k, ls in enumerate(x.serial for x in at.ligands):
                    if ls not in serials:
                        at.is_root = True
                        rooted = True
                        break
                if rooted:
                    break
            assert rooted, 'Could not find a root atom for a glycan -- this is a bug'
            out.append((np.array([x.serial - 1 for x in g], dtype=int),) + self._rotatable_bonds(g))
        return out

    @staticmethod
    def _rotatable_bonds(g):
        """From a rooted pendant graph (one or more atoms flagged ``is_root`` where the group is
        anchored to the rest of the structure), return ``(bonds, movers)``: each bridging bond that
        is not a bond-to-H and not a peptide bond, cut into (bond endpoints, unrooted-side mover
        atoms).  A bond whose *both* sides carry a root (e.g. an internal backbone bond of a loop
        anchored at both ends) is not rotatable and is skipped."""
        bonds, movers = [], []
        for ai, aj in list(nx.bridges(g)):
            if (ai.isH() or aj.isH()) or ai.is_pep(aj):
                continue
            g.remove_edge(ai, aj)
            S = [g.subgraph(c).copy() for c in nx.connected_components(g)]
            assert len(S) == 2, f'Bond {ai.serial-1}-{aj.serial-1} is not a bridge'
            rooted = [any(getattr(x, 'is_root', False) for x in sg) for sg in S]
            if False in rooted:                       # exactly one unrooted side -> rotatable
                unrooted = S[rooted.index(False)]
                unrooted.remove_nodes_from([n for n in (ai, aj) if n in unrooted])
                nodes = list(unrooted.nodes)
                if len(nodes) > 1 or (len(nodes) == 1 and not nodes[0].isH()):
                    bonds.append((ai.serial - 1, aj.serial - 1))
                    movers.append(np.array([x.serial - 1 for x in nodes], dtype=int))
            g.add_edge(ai, aj)
        return bonds, movers

    def _run_pendant_declash(self, state, pendants, cycles, clashdist):
        """Load the current psf/pdb, run the numpy pendant declasher over each ``(indices, bonds,
        movers)`` entry (sharing one seeded RNG), and write ``{basename}.pdb``.

        The clash environment is the *static* set of heavy atoms outside every pendant in the batch
        (e.g. the protein, for glycans).  Excluding the other movable pendants keeps the sequential
        greedy from chasing them into each other; the destabilizing overlaps that matter for the
        downstream minimization are pendant-vs-structure, and inter-pendant contacts relax out.
        """
        from ..molecule.coordmanip import CoordManipulator
        from ..psfutil.declash import declash_pendant
        cm = CoordManipulator(state.psf.name, state.pdb.name)
        coords = cm.coords
        heavy = cm._mass > 1.1                        # PSF masses: hydrogens ~1.008
        batch = np.zeros(len(coords), dtype=bool)
        for indices, _bonds, _movers in pendants:
            batch[indices] = True
        env_heavy = np.nonzero(heavy & ~batch)[0]
        rng = np.random.default_rng(self._DECLASH_SEED)
        for indices, bonds, movers in pendants:
            pend = np.zeros(len(coords), dtype=bool)
            pend[indices] = True
            pendant_heavy = np.nonzero(heavy & pend)[0]
            declash_pendant(coords, pendant_heavy, env_heavy, bonds, movers, cycles, clashdist, rng)
        cm.coords = coords
        cm.write_pdb(f'{self.basename}.pdb')

    #: fixed RNG seed so declashed builds are reproducible (VMD's rand() was unseeded)
    _DECLASH_SEED = 90125

    def ingest_molecules(self):
        """
        Ingests the base molecule from the specifications provided in the task.
        This molecule is initialized from either base coordinates fetched in a prior fetch task,
        or a state defined by a previous continuation task.
        It also handles any graft sources specified in the sequence modifications and
        activates the biological assembly of the base molecule.
        """
        specs = self.specs
        self.source_specs = specs['source']
        assert not 'id' in self.source_specs, f'Version 2.0+ of Pestifer does not support "id" in source specs.  Please begin your task list with either a `fetch` or `continuation` task.'
        this_source = {}
        base_coordinates: Path = self.get_current_artifact_path('base_coordinates')
        if not base_coordinates:
            state: StateArtifacts = self.get_current_artifact('state')
            if not (state.pdb and state.psf):
                raise PestiferBuildError(f'No base_coordinates artifact found, and no prebuilt PDB/PSF/XSC files found in the pipeline context. Cannot ingest base molecule.')
            this_source['prebuilt'] = {
                'pdb': state.pdb.name,
                'psf': state.psf.name,
                'xsc': state.xsc.name if state.xsc else None
            }
            this_source['file_format'] = 'PDB'
        else:
            basename, ext = os.path.splitext(base_coordinates)
            this_source['source_id'] = basename
            this_source['source_db'] = 'local'
            if ext.lower() == '.pdb':
                this_source['file_format'] = 'PDB'
            elif ext.lower() == '.cif':
                this_source['file_format'] = 'mmCIF'
            else:
                raise PestiferBuildError(f'Unknown file format {ext} for base_coordinates artifact {base_coordinates.name}')
        self.source_specs.update(this_source)
        logger.debug(f'User-input modspecs:')
        my_logger(self.specs["mods"], logger.debug, depth=1)
        self.objmanager = ObjManager(self.specs['mods'])
        seqmods = self.objmanager.get('seq', {})
        logger.debug(f'Ingesting seqmods:')
        my_logger(seqmods, logger.debug, depth=1)
        if 'grafts' in seqmods:
            # logger.debug(f'looking for graft sources to ingest')
            Grafts: GraftList = seqmods['grafts']
            graft_artifacts = []
            for g in Grafts.data:
                if not g.source_pdbid in self.molecules:
                    logger.debug(f'Ingesting graft source {g.source_pdbid}')
                    this_source = {
                        'source_db': 'rcsb',
                        'source_id': g.source_pdbid,
                        'file_format': 'PDB',
                        'sequence': {'skip_glycan_renaming': True}
                    }
                    # the Molecule call below will fetch coordinates for graft pdbs
                    self.molecules[g.source_pdbid] = Molecule(source=this_source)
                    graft_artifacts.append(PDBFileArtifact(g.source_pdbid))
                g.activate(deepcopy(self.molecules[g.source_pdbid]))
            if len(graft_artifacts) > 0:
                self.register(graft_artifacts, key='graft_sources', artifact_type=PDBFileArtifactList)
        if 'Cfusions' in seqmods:
            # A Cfusion's donor coordinates are loaded directly in VMD by the psfgen scripter
            # (`mol new <sourcefile>`), so we only need the raw PDB present -- not a parsed
            # Molecule (parsing it would also raise spurious inter-residue link warnings for a
            # structure we only mine coordinates from).  If the sourcefile isn't already local,
            # treat its basename as an RCSB PDB ID and fetch the file (parity with grafts, which
            # name a source PDB ID), so a Cfusion can name a PDB ID and pestifer provisions it.
            Cfusions: CfusionList = seqmods['Cfusions']
            cfusion_artifacts = []
            for c in Cfusions.data:
                pdbid = os.path.splitext(os.path.basename(c.sourcefile))[0]
                target = f'{pdbid}.pdb'
                if not os.path.exists(c.sourcefile):
                    if not os.path.exists(target):
                        logger.debug(f'Fetching Cfusion donor {pdbid}')
                        PDBParser(source_db='rcsb', source_id=pdbid).fetch()
                        cfusion_artifacts.append(PDBFileArtifact(pdbid))
                    c.sourcefile = target
            if len(cfusion_artifacts) > 0:
                self.register(cfusion_artifacts, key='cfusion_sources', artifact_type=PDBFileArtifactList)
        self.chainIDmanager = ChainIDManager(
            format = self.source_specs['file_format'],
            transform_reserves = self.source_specs.get('transform_reserves', {}),
            remap = self.source_specs.get('remap_chainIDs', {}))
        # if 'vmdatomselections' in self.source_specs['exclude']:
        #     self.apply_vmdexclusions()
        self.base_molecule = Molecule(source = self.source_specs,
                                      objmanager = self.objmanager,
                                      chainIDmanager = self.chainIDmanager).activate_biological_assembly(self.source_specs['biological_assembly'])
        # register self.base_molecule in the pipeline context
        self.register(self.base_molecule, key='base_molecule')
        for molid, molecule in self.molecules.items():
            logger.debug(f'Molecule "{molid}": {molecule.num_atoms()} atoms in {molecule.num_residues()} residues; {molecule.num_segments()} segments.')
        if 'cif' in self.source_specs.get('file_format', '').lower():
            if self.source_specs.get('cif_residue_map_file', ''):
                write_residue_map(self.base_molecule.asymmetric_unit.residues.cif_residue_map(), self.source_specs['cif_residue_map_file'])
                self.register(Path(self.source_specs['cif_residue_map_file']), key='cif_residue_map', artifact_type=CSVDataFileArtifact)
            logger.debug(f'Base molecule has {self.base_molecule.num_atoms()} atoms in {self.base_molecule.num_residues()} residues; {self.base_molecule.num_segments()} segments.')

    def update_molecule(self):
        """
        Updates all segments of the base molecule based on the 
        current coordinate file.  All ssbonds and links are 
        carried forward. No biological assembly beyond the apparent
        asymmetric unit is assumed. This should permit generation
        of a new psfgen script based on this new molecule to 
        recreate it or modify it.
        """
        # get the key of the base_molecule
        logger.debug(f'{self.taskname} has {len(self.molecules)} entries in its molecules dict')
        base_key = 'base'
        for k, v in self.molecules.items():
            if v == self.base_molecule:
                base_key = k
        # assert base_key!='UNSET',f'Cannot update a non-existent base molecule'
        # get the psf, pdb, xsc from the pipeline context
        state: StateArtifacts = self.get_current_artifact('state')
        pdb: Path = state.pdb
        psf: Path = state.psf
        xsc: Path | None = state.xsc
        if not (pdb and psf):
            raise PestiferBuildError(f'No base_coordinates artifact found, and no prebuilt PDB/PSF/XSC files found in the pipeline context. Cannot ingest base molecule.')
        source = {}
        source['prebuilt'] = {
            'pdb': pdb.name,
            'psf': psf.name,
            'xsc': xsc.name if xsc else None
        }
        if hasattr(self, 'chainIDmanager') and hasattr(self, 'objmanager'):
            updated_molecule = Molecule(source=source, chainIDmanager=self.chainIDmanager, objmanager=self.objmanager).activate_biological_assembly(0)
        else:
            updated_molecule = Molecule(source=source).activate_biological_assembly(0)

        self.molecules[base_key] = updated_molecule
        self.base_molecule = updated_molecule
        self.register(self.base_molecule, key='base_molecule')
