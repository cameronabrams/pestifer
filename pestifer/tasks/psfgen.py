# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`PsfgenTask` class for handling invocations of psfgen which create a molecule from a base PDB/mmCIF file.

Usage is described in the :ref:`subs_runtasks_psfgen` documentation.

"""
import logging
import networkx as nx
import os

from copy import deepcopy
from pathlib import Path
from typing import ClassVar


from .basetask import VMDTask
from ..molecule.chainidmanager import ChainIDManager
from ..molecule.molecule import Molecule
from ..core.objmanager import ObjManager
from ..core.artifacts import *
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
    This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class, and as such it has two parameters.

    Parameters
    ----------
    config_specs : dict
        Configuration specifications for the task.
    """

    _yaml_header: ClassVar[str] = 'psfgen'
    """
    YAML header for the PsfgenTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """

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
                if len(objlist) > 0:
                    self.next_basename(objtype)
                    vm: VMDScripter = self.scripters['vmd']
                    packages = []
                    if objtype == 'crotations':
                        packages.append('PestiferCRot')
                    vm.newscript(self.basename, packages=packages)
                    state: StateArtifacts = self.get_current_artifact('state')
                    vm.load_psf_pdb(state.psf.name, state.pdb.name, new_molid_varname='mCM')
                    match objtype:
                        case 'crotations':
                            for transform in self.base_molecule.active_biological_assembly.transforms.data:
                                vm.write_crots(objlist, chainIDmap=transform.chainIDmap)
                        case 'orient':
                            vm.write_orients(objlist)
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
                raise ValueError(f'No topology file found for residue name {resname}')
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
        return 0
        
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
        Declash loops in the base molecule using the custom ``PestiferDeclash`` TcL package.
        This method generates a VMD script to identify and declash loops in the molecular structure.
        It uses the ``PestiferDeclash`` package to perform the declashing operation, which involves identifying
        and modifying the coordinates of atoms in loop regions.
        """
        specs_to_declashers = specs#['source']['sequence']['declash']
        if segtype == 'protein':
            self.declash_protein_loops(specs_to_declashers)
        elif segtype == 'nucleicacid':
            self.declash_na_loops(specs_to_declashers)
        elif segtype == 'glycan':
            self.declash_glycans(specs_to_declashers)
        else:
            raise ValueError(f'Unknown segment type {segtype} for declashing')

    def declash_protein_loops(self, specs: dict):
        """
        Declash loops in the base molecule using the custom ``PestiferDeclash`` TcL package.
        This method generates a VMD script to identify and declash loops in the molecular structure.
        It uses the ``PestiferDeclash`` package to perform the declashing operation, which involves identifying
        and modifying the coordinates of atoms in loop regions.
        """
        mol = self.base_molecule
        cycles = specs['maxcycles']
        self.next_basename('declash-loops')
        state: StateArtifacts = self.get_current_artifact('state')
        psf = state.psf.name
        pdb = state.pdb.name
        vt: VMDScripter = self.scripters['vmd']
        vt.newscript(self.basename, packages=['PestiferDeclash'])
        vt.load_psf_pdb(psf, pdb, new_molid_varname='mLL')
        vt.write_protein_loop_lines(mol, cycles=cycles, include_c_termini=specs['include_C_termini'])
        vt.write_pdb(self.basename, 'mLL')
        vt.writescript()
        vt.runscript()
        self.register(dict(
            pdb=PDBFileArtifact(self.basename), 
            psf=state.psf), key='state', artifact_type=StateArtifacts)
        self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
        self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)

    def declash_na_loops(self,specs):
        """
        Declash nucleic acid loops in the base molecule using the custom ``PestiferDeclash`` TcL package.
        This method generates a VMD script to identify and declash nucleic acid loops in the molecular structure.
        It uses the ``PestiferDeclash`` package to perform the declashing operation, which involves identifying
        and modifying the coordinates of atoms in nucleic acid loop regions.
        """
        cycles = specs['maxcycles']
        clashdist = specs['clashdist']
        self.next_basename('declash-na-loops')
        vt: VMDScripter = self.scripters['vmd']
        state: StateArtifacts = self.get_current_artifact('state')
        psf = state.psf.name
        pdb = state.pdb.name
        outpdb = f'{self.basename}.pdb'
        vt.newscript(self.basename, packages=['PestiferDeclash'])
        vt.addline(f'mol new {psf}')
        vt.addline(f'mol addfile {pdb} waitfor all')
        vt.addline(f'set a [atomselect top all]')
        vt.addline(f'set molid [molinfo top get id]')
        nna = self._write_na_loops(vt)
        vt.addline(f'set nna {nna}')
        vt.addline(f'vmdcon -info "Declashing $nna nucleic acid loops; clashdist {clashdist}; maxcycles {cycles}"')
        vt.addline(r'for {set i 0} {$i<$nna} {incr i} {')
        vt.addline(f'   declash_pendant $molid $na_idx($i) $rbonds($i) $movers($i) {cycles} {clashdist}')
        vt.addline(r'}')
        vt.addline(f'$a writepdb {outpdb}')
        vt.writescript()
        logger.debug(f'Declashing {nna} nucleic acid loops')
        vt.runscript(progress_title='declash-nucleic-acid-loops')
        self.register(dict(
            pdb=PDBFileArtifact(self.basename), 
            psf=state.psf), key='state', artifact_type=StateArtifacts)
        self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
        self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)

    def _write_na_loops(self, vt: VMDScripter, **options):
        mol = self.base_molecule
        au = mol.asymmetric_unit
        state: StateArtifacts = self.get_current_artifact('state')
        psf = state.psf.name
        logger.debug(f'ingesting {psf}')
        struct = PSFContents(psf, parse_topology=['bonds'])
        logger.debug(f'PSF has {len(struct.atoms)} atoms; extracting nucleic acid atoms...')
        na_atoms: PSFAtomList = struct.atoms.get(lambda x: x.segtype == 'nucleicacid')
        logger.debug(f'PSF has {len(na_atoms)} nucleic acid atoms')
        # my_rep = list(set([(x.chainID, x.resid) for x in na_atoms.data]))
        # my_rep.sort(key=lambda x: (x[0], x[1]))
        logger.debug(f'Getting loops from {len(na_atoms)} nucleic acid atoms in PSF file {psf}')
        # logger.debug(f'{my_rep}')
        min_length = self.min_loop_length
        include_c_termini = options.get('include_c_termini', False)
        i = 0
        SL = [S for S in au.segments.data if S.segtype == 'nucleicacid']
        logger.debug(f'Asymmetric unit has {len(SL)} nucleic acid segments')
        for S in SL:
            asymm_segname = S.segname
            logger.debug(f'Processing segment {asymm_segname} for nucleic acid loops')
            n_subsegs = len(S.subsegments)
            for b in S.subsegments.data:
                lr_resid = S.residues[b.bounds[0]].resid
                rr_resid = S.residues[b.bounds[1]].resid
                logger.debug(f'Processing subsegment {b.state} for segname {asymm_segname} with bounds {lr_resid}-{rr_resid}')
                is_c_terminus = (S.subsegments.index(b) == (n_subsegs - 1))
                is_processible = b.state == 'MISSING' and b.num_items() >= min_length
                if is_processible and (not include_c_termini) and is_c_terminus:
                    logger.debug(f'A.U. C-terminal loop {b.state} declashing is skipped')
                    is_processible = False
                if is_processible:
                    logger.debug(f'Processing loop {b.state} {b.bounds} for segname {asymm_segname}')
                    loop_atoms = PSFAtomList([x for x in na_atoms.data if x.segname == asymm_segname and x.resid >= lr_resid and x.resid <= rr_resid])
                    logger.debug(f'Loop {b.state} has {len(loop_atoms)} atoms from PSFAtomList')
                    na_graph = loop_atoms.graph()
                    logger.debug(f'{na_graph}')
                    G = [na_graph.subgraph(c).copy() for c in nx.connected_components(na_graph)]
                    assert len(G) == 1, f'NA loop {b.state} has more than one connected component'
                    logger.debug(f'Loop {b.state} has {len(loop_atoms)} atoms')
                    g = G[0]
                    serials = [x.serial for x in g]
                    for at in g:
                        lig_ser = [x.serial for x in at.ligands]
                        for k, ls in enumerate(lig_ser):
                            if not ls in serials:
                                at.is_root = True
                                rp = at.ligands[k]
                                logger.debug(f'-> Atom {str(at)} is the root, bound to atom {str(rp)}')
                    indices = ' '.join([str(x.serial-1) for x in g])
                    vt.addline(f'set na_idx({i}) [list {indices}]')
                    vt.addline(f'set rbonds({i}) [list]')
                    vt.addline(f'set movers({i}) [list]')
                    for bond in nx.bridges(g):
                        ai, aj = bond
                        if not (ai.isH() or aj.isH()) and not ai.is_pep(aj):
                            g.remove_edge(ai, aj)
                            CC = [g.subgraph(c).copy() for c in nx.connected_components(g)]
                            assert len(CC) == 2, f'Bond {ai.serial-1}-{aj.serial-1} when cut makes more than 2 components'
                            for sg in CC:
                                is_root = any([hasattr(x, 'is_root') for x in sg])
                                if not is_root:
                                    if ai in sg:
                                        sg.remove_node(ai)
                                    if aj in sg:
                                        sg.remove_node(aj)
                                    if len(sg) > 1 or (len(sg) == 1 and not [x for x in sg.nodes][0].isH()):
                                        mover_serials = [x.serial for x in sg]
                                        mover_indices = " ".join([str(x-1) for x in mover_serials])
                                        logger.debug(f'{str(ai)}--{str(aj)} is a rotatable bridging bond')
                                        vt.addline(f'lappend rbonds({i}) [list {ai.serial-1} {aj.serial-1}]')
                                        logger.debug(f'  -> movers: {" ".join([str(x) for x in sg])}')
                                        vt.addline(f'lappend movers({i}) [list {mover_indices}]')
                            g.add_edge(ai, aj)
                    i += 1
        return i

    def declash_glycans(self, specs):
        """
        Declash glycans in the base molecule using the custom ``PestiferDeclash`` TcL package.
        This method generates a VMD script to identify and declash glycans in the molecular structure.
        It uses the ``PestiferDeclash`` package to perform the declashing operation,
        which involves identifying and modifying the coordinates of atoms in glycan regions.
        """
        mol: Molecule = self.base_molecule
        cycles: int = specs['maxcycles']
        clashdist: float = specs['clashdist']
        if not mol.nglycans() or not cycles:
            logger.debug(f'Glycan declashing is intentionally not done.')
            return
        self.next_basename('declash-glycans')
        outpdb = f'{self.basename}.pdb'
        state: StateArtifacts = self.get_current_artifact('state')
        psf = state.psf.name
        pdb = state.pdb.name
        vt: VMDScripter = self.scripters['vmd']
        vt.newscript(self.basename, packages=['PestiferDeclash'])
        vt.addline(f'mol new {psf}')
        vt.addline(f'mol addfile {pdb} waitfor all')
        vt.addline(f'set a [atomselect top all]')
        vt.addline(f'set molid [molinfo top get id]')
        nglycan = self._write_glycans(vt)
        vt.addline(f'vmdcon -info "Declashing $nglycan glycans; clashdist {clashdist}; maxcycles {cycles}"')
        vt.addline(r'for {set i 0} {$i<$nglycan} {incr i} {')
        vt.addline(f'   declash_pendant $molid $glycan_idx($i) $rbonds($i) $movers($i) {cycles} {clashdist}')
        vt.addline(r'}')
        vt.addline(f'$a writepdb {outpdb}')
        vt.writescript()
        logger.debug(f'Declashing {nglycan} glycans')
        vt.runscript(progress_title='declash-glycans')
        self.register(dict(
            pdb=PDBFileArtifact(self.basename), 
            psf=state.psf), key='state', artifact_type=StateArtifacts)
        self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
        self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)

    def _write_glycans(self, fw: VMDScripter):
        state: StateArtifacts = self.get_current_artifact('state')
        psf = state.psf.name
        logger.debug(f'ingesting {psf}')
        struct = PSFContents(psf, parse_topology=['bonds'])
        logger.debug(f'Making graph structure of glycan atoms...')
        glycanatoms = struct.atoms.get(lambda x: x.segtype == 'glycan')
        logger.debug(f'{len(glycanatoms)} total glycan atoms')
        glycangraph = glycanatoms.graph()
        G = [glycangraph.subgraph(c).copy() for c in nx.connected_components(glycangraph)]
        logger.debug(f'Preparing declash input for {len(G)} glycans')
        fw.addline(f'set nglycans {len(G)}')
        for i, g in enumerate(G):
            logger.debug(f'Glycan {i} has {len(g)} atoms')
            serials = [x.serial for x in g]
            for j, at in enumerate(g):
                lig_ser = [x.serial for x in at.ligands]
                # logger.debug(f'Atom {at.serial} ({j}) has {len(lig_ser)} ligands')
                for k, ls in enumerate(lig_ser):
                    if not ls in serials:
                        at.is_root = True
                        rp = at.ligands[k]
                        # logger.debug(f'-> Atom {at.serial} ({j}) is the root, bound to atom {rp.serial}')
            indices = ' '.join([str(x.serial-1) for x in g])
            # logger.debug(f'indices {indices}')
            fw.comment(f'Glycan {i}:')
            fw.addline(f'set glycan_idx({i}) [list {indices}]')
            fw.addline(f'set rbonds({i}) [list]')
            fw.addline(f'set movers({i}) [list]')
            for bond in nx.bridges(g):
                ai, aj = bond
                if not (ai.isH() or aj.isH()) and not ai.is_pep(aj):
                    g.remove_edge(ai, aj)
                    S = [g.subgraph(c).copy() for c in nx.connected_components(g)]
                    assert len(S) == 2, f'Bond {ai.serial-1}-{aj.serial-1} when cut makes more than 2 components'
                    for sg in S:
                        is_root = any([hasattr(x, 'is_root') for x in sg])
                        if not is_root:
                            if ai in sg:
                                sg.remove_node(ai)
                            if aj in sg:
                                sg.remove_node(aj)
                            if len(sg) > 1 or (len(sg) == 1 and not [x for x in sg.nodes][0].isH()):
                                mover_serials = [x.serial for x in sg]
                                mover_indices = " ".join([str(x-1) for x in mover_serials])
                                logger.debug(f'{str(ai)}--{str(aj)} is a rotatable bridging bond')
                                fw.addline(f'lappend rbonds({i}) [list {ai.serial-1} {aj.serial-1}]')
                                logger.debug(f'  -> movers: {" ".join([str(x) for x in sg])}')
                                fw.addline(f'lappend movers({i}) [list {mover_indices}]')
                    g.add_edge(ai, aj)
        return len(G)

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
                raise RuntimeError(f'No base_coordinates artifact found, and no prebuilt PDB/PSF/XSC files found in the pipeline context. Cannot ingest base molecule.')
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
                raise RuntimeError(f'Unknown file format {ext} for base_coordinates artifact {base_coordinates.name}')
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
                        'file_format': 'PDB'
                    }
                    # the Molecule call below will fetch coordinates for graft pdbs
                    self.molecules[g.source_pdbid] = Molecule(source=this_source)
                    graft_artifacts.append(PDBFileArtifact(g.source_pdbid))
                g.activate(deepcopy(self.molecules[g.source_pdbid]))
            if len(graft_artifacts) > 0:
                self.register(graft_artifacts, key='graft_sources', artifact_type=PDBFileArtifactList)
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
            raise RuntimeError(f'No base_coordinates artifact found, and no prebuilt PDB/PSF/XSC files found in the pipeline context. Cannot ingest base molecule.')
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
