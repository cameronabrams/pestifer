# Author: Cameron F. Abrams, <cfa2@drexel.edu>
"""
Definition of the :class:`MakeMembraneSystemTask` class for handling embedding proteins into bilayers.

Usage is described in the :ref:`subs_runtasks_make_membrane_system` documentation.
"""

import glob
import logging
import numpy as np
from pathlib import Path

from .basetask import BaseTask
# from .terminate import TerminateTask

from ..charmmff.charmmffcontent import CHARMMFFContent

from ..core.artifacts import *
from ..util.stringthings import __pestifer_version__

from ..molecule.bilayer import Bilayer, specstrings_builddict

from ..psfutil.psfcontents import get_toppar_from_psf

from ..scripters import PsfgenScripter, PackmolScripter, VMDScripter

from ..util.util import cell_to_xsc,cell_from_xsc, protect_str_arg
from ..util.units import _UNITS_

sA2_ = _UNITS_['SQUARE-ANGSTROMS']

logger = logging.getLogger(__name__)

class MakeMembraneSystemTask(BaseTask):
    """ 
    A class for handling embedding proteins into bilayers
    """
    _yaml_header = 'make_membrane_system'
    """
    YAML header for the MakeMembraneSystemTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """

    def provision(self, packet: dict):
        logger.debug(f'Provisioning MakeMembraneSystemTask with packet:')
        my_logger(packet, logger.debug)
        super().provision(packet)
        self.patchA: Bilayer = None
        self.patchB: Bilayer = None
        self.patch: Bilayer = None
        self.bilayer_specs: dict = self.specs.get('bilayer', {})
        self.embed_specs: dict = self.specs.get('embed', {})
        self.using_prebuilt_bilayer: bool = False
        self.progress = self.provisions.get('progress-flag', True)
        self.charmmff_content: CHARMMFFContent = self.resource_manager.charmmff_content
    
        if 'prebuilt' in self.bilayer_specs and 'pdb' in self.bilayer_specs['prebuilt']:
            logger.debug('Using prebuilt bilayer')
            self.using_prebuilt_bilayer = True
            self.quilt: Bilayer = Bilayer()
            quilt_state = dict(pdb=self.bilayer_specs['prebuilt']['pdb'], psf=self.bilayer_specs['prebuilt']['psf'], xsc=self.bilayer_specs['prebuilt']['xsc'])
            self.register(quilt_state, key='quilt_state', artifact_type=StateArtifacts)
            self.quilt.box, self.quilt.origin = cell_from_xsc(quilt_state.xsc.path)
            self.quilt.area = self.quilt.box[0][0] * self.quilt.box[1][1]
            additional_topologies = get_toppar_from_psf(quilt_state.psf.name)
            # these will be registered as artifacts when psfgen executes
            self.quilt.addl_streamfiles = additional_topologies
        else:
            self.initialize()

    def initialize(self):
        """
        Initialize the MakeMembraneSystemTask by building the bilayer patch and quilt.
        This method sets up the bilayer patch based on the specifications provided in the configuration.
        If the bilayer is asymmetric, it builds two symmetric patches.
        It also sets up the quilt from the bilayer patch, which will be used for embedding proteins. 
        If a prebuilt bilayer is specified, it uses that instead of building a new one.
        """
        lipid_specstring = self.bilayer_specs.get('lipids', '')
        ratio_specstring = self.bilayer_specs.get('mole_fractions', '')
        conformers_specstring = self.bilayer_specs.get('conformers', '')
        solvent_specstring = self.bilayer_specs.get('solvents', 'TIP3')
        solvent_ratio_specstring = self.bilayer_specs.get('solvent_mole_fractions', '1.0')
        solvent_to_lipid_ratio = self.bilayer_specs.get('solvent_to_lipid_ratio', 32.0)
        patch_nlipids = self.bilayer_specs.get('patch_nlipids', dict(upper=100, lower=100))
        cation_name = self.bilayer_specs.get('cation', 'POT')
        anion_name = self.bilayer_specs.get('anion', 'CLA')
        neutralizing_salt = [cation_name, anion_name]
        salt_con = self.bilayer_specs.get('salt_con', 0.0)  # Molar concentration
        composition_dict = self.bilayer_specs.get('composition', {})

        if not composition_dict['upper_leaflet'] or not composition_dict['lower_leaflet']:
            logger.debug('No upper or lower leaflet specified in composition; building from memgen-format specstrings')
            composition_dict = specstrings_builddict(lipid_specstring,
                                                     ratio_specstring,
                                                     conformers_specstring,
                                                     solvent_specstring,
                                                     solvent_ratio_specstring)
        logger.debug(f'Naive main composition dict:')
        my_logger(composition_dict, logger.debug)
        self.patch = Bilayer(composition_dict,
                             neutralizing_salt = neutralizing_salt,
                             salt_concentration = salt_con,
                             solvent_specstring = solvent_specstring,
                             solvent_ratio_specstring = solvent_ratio_specstring,
                             solvent_to_key_lipid_ratio = solvent_to_lipid_ratio,
                             leaflet_nlipids = patch_nlipids,
                             charmmffcontent = self.charmmff_content)
        for spdb in self.patch.register_species_pdbs:
            self.register(spdb, key='species_pdb_for_packmol', artifact_type=PDBFileArtifact)
        logger.debug(f'Mature main composition dict:')
        my_logger(composition_dict, logger.debug)
        if self.patch.asymmetric:
            logger.debug(f'Requested patch is asymmetric; generating two symmetric patches')
            logger.debug(f'Symmetrizing bilayer to upper leaflet')
            composition_dict['lower_leaflet_saved'] = composition_dict['lower_leaflet']
            composition_dict['lower_chamber_saved'] = composition_dict['lower_chamber']
            composition_dict['lower_leaflet'] = composition_dict['upper_leaflet']
            composition_dict['lower_chamber'] = composition_dict['upper_chamber']
            self.patchA = Bilayer(composition_dict,
                                  neutralizing_salt = neutralizing_salt,
                                  salt_concentration = salt_con,
                                  solvent_specstring = solvent_specstring,
                                  solvent_ratio_specstring = solvent_ratio_specstring,
                                  solvent_to_key_lipid_ratio = solvent_to_lipid_ratio,
                                  leaflet_nlipids = patch_nlipids,
                                  charmmffcontent = self.charmmff_content)
            logger.debug(f'Symmetrizing bilayer to lower leaflet')
            composition_dict['upper_leaflet_saved'] = composition_dict['upper_leaflet']
            composition_dict['upper_chamber_saved'] = composition_dict['upper_chamber']
            composition_dict['lower_leaflet'] = composition_dict['lower_leaflet_saved']
            composition_dict['lower_chamber'] = composition_dict['lower_chamber_saved']
            composition_dict['upper_leaflet'] = composition_dict['lower_leaflet']
            composition_dict['upper_chamber'] = composition_dict['lower_chamber']
            self.patchB = Bilayer(composition_dict,
                                  neutralizing_salt = neutralizing_salt,
                                  solvent_specstring = solvent_specstring,
                                  solvent_ratio_specstring = solvent_ratio_specstring,
                                  solvent_to_key_lipid_ratio = solvent_to_lipid_ratio,
                                  leaflet_nlipids = patch_nlipids,
                                  charmmffcontent = self.charmmff_content)
            composition_dict['upper_leaflet'] = composition_dict['upper_leaflet_saved']
            composition_dict['upper_chamber'] = composition_dict['upper_chamber_saved']
            self.patch = None

    def do(self) -> int:
        """
        Execute the MakeMembraneSystemTask.
        """
        # as part of a list of tasks, this task expects to be fed a protein system to embed
        protein_state: StateArtifacts = self.get_current_artifact('state')
        if protein_state is None or protein_state.psf is None or protein_state.pdb is None:
            self.embedding = False
        else:
            self.embedding = True
            logger.debug(f'Using psf {protein_state.psf.name} and pdb {protein_state.pdb.name} as inputs')

        if self.embedding:
            self.orient_protein()
        if not self.using_prebuilt_bilayer:
            self.build_patch()
            self.make_quilt_from_patch()
        if self.embedding:
            self.embed_protein()
        else:
            self.pipeline.rekey('quilt_state', 'state')
        return 0

    def orient_protein(self):
        """
        Orient the protein so that its principal axis is aligned with the z-axis (membrane normal).
        This method uses the Orient package in VMD to perform the orientation based on specified atom selections.
        The oriented protein structure is then saved as a new PDB file and registered as an artifact.
        """
        if self.embed_specs.get('no_orient', True):
            logger.debug('Skipping orientation of protein as per embed specs')
            return
        
        z_head_group = self.embed_specs.get('z_head_group', None)
        z_tail_group = self.embed_specs.get('z_tail_group', None)
        if z_head_group is None or z_tail_group is None:
            logger.debug('No z_head_group or z_tail_group specified; skipping orientation of protein')
            return
        
        self.next_basename('orient')
        logger.debug(f'Orienting protein')
        vm: VMDScripter = self.scripters['vmd']
        vm.newscript(self.basename)
        vm.usescript('bilayer_orient')
        vm.writescript(self.basename)
        self.register(self.basename, key='tcl', artifact_type=VMDScriptArtifact)
        protein_state: StateArtifacts = self.get_current_artifact('state')
        protein_psf: str = protein_state.psf.name
        protein_pdb: str = protein_state.pdb.name
        result = vm.runscript(psf=protein_psf,
                              pdb=protein_pdb,
                              z_head_group=protect_str_arg(z_head_group),
                              z_tail_group=protect_str_arg(z_tail_group),
                              o=self.basename)
        if result != 0:
            raise RuntimeError(f'vmd failed with result {result} for {self.basename}')
        self.register(self.basename, key='log', artifact_type=VMDLogFileArtifact)
        self.register(dict(
            psf=protein_psf, 
            pdb=PDBFileArtifact(self.basename, pytestable=True)), 
            key='state', artifact_type=StateArtifacts)

    def build_patch(self):
        """
        Build the bilayer patch or patches based on the specifications provided in the configuration.
        This method retrieves the bilayer specifications, including solution conditions, rotation parameters,
        and other relevant settings.
        It then constructs the patch or patches, packs them using Packmol, generates the PSF file, and then does a short series of equilibration MD simulations.
        """
        logger.debug(f'Bilayer specs:')
        my_logger(self.bilayer_specs, logger.debug)
        solution_gcc: float = self.bilayer_specs.get('solution_gcc',1.0)
        rotation_pm: float = self.bilayer_specs.get('rotation_pm',10.)
        half_mid_zgap: float = self.bilayer_specs.get('half_mid_zgap',1.0)
        SAPL: float = self.bilayer_specs.get('SAPL',75.0)
        seed: int = self.bilayer_specs.get('seed',27021972)
        tolerance: float = self.bilayer_specs.get('tolerance',2.0)
        xy_aspect_ratio: float = self.bilayer_specs.get('xy_aspect_ratio',1.0)
        nloop: int = self.bilayer_specs.get('nloop',100)
        nloop_all: int = self.bilayer_specs.get('nloop_all',100)
        relaxation_protocols: dict = self.bilayer_specs.get('relaxation_protocols',{})
        relaxation_protocol: dict = relaxation_protocols.get('patch',{})
        logger.debug('Relaxation protocols:')
        my_logger(relaxation_protocols, logger.debug)
        # we now build the patch, or if asymmetric, two patches
        for patch, specbyte in zip([self.patch, self.patchA, self.patchB], ['', 'A', 'B']):
            if patch is None:
                continue
            patch.spec_out(SAPL=SAPL, xy_aspect_ratio=xy_aspect_ratio,
                            rotation_pm=rotation_pm, solution_gcc=solution_gcc,
                            half_mid_zgap=half_mid_zgap)
            self.pack_patch(patch, patch_name=f'patch{specbyte}', seed=seed,
                            tolerance=tolerance,
                            nloop_all=nloop_all,
                            half_mid_zgap=half_mid_zgap,
                            rotation_pm=rotation_pm,
                            nloop=nloop)
            self.do_psfgen(patch, bilayer_name=f'patch{specbyte}')
            self.equilibrate_bilayer(patch, bilayer_name=f'patch{specbyte}', relaxation_protocol=relaxation_protocol)

    def register_tops_streams_from_psfgen(self, filelist):
        self.register([CharmmffTopFileArtifact(x) for x in filelist if x.endswith('rtf')], key='charmmff_topfiles', artifact_type=CharmmffTopFileArtifacts)
        self.register([CharmmffStreamFileArtifact(x) for x in filelist if x.endswith('str')], key='charmmff_streamfiles', artifact_type=CharmmffStreamFileArtifacts)

    def register_tmpfiles_from_tcl(self, logname: str):
        tmp_files = {}
        with open(logname, 'r') as f:
            lines = f.readlines()
            for line in lines:
                if 'register_as_artifact' in line:
                    artifactname = line.split()[-1]
                    ext = os.path.splitext(artifactname)[1].replace('.','')
                    if not ext in tmp_files:
                        tmp_files[ext] = []
                    tmp_files[ext].append(artifactname)
        logger.debug(f'tmp_files: {tmp_files}')
        if 'pdb' in tmp_files:
            self.register([PDBFileArtifact(x) for x in tmp_files['pdb']], key='tmp_patch_pdbs', artifact_type=PDBFileArtifactList)
        if 'psf' in tmp_files:
            self.register([PSFFileArtifact(x) for x in tmp_files['psf']], key='tmp_patch_psfs', artifact_type=PSFFileArtifactList)
        if 'coor' in tmp_files:
            self.register([NAMDCoorFileArtifact(x) for x in tmp_files['coor']], key='tmp_patch_coors', artifact_type=NAMDCoorFileArtifactList)
        if 'log' in tmp_files:
            self.register([NAMDLogFileArtifact(x) for x in tmp_files['log']], key='tmp_patch_logs', artifact_type=NAMDLogFileArtifactList)

    def do_psfgen(self, patch: Bilayer, bilayer_name: str):
        """
        Perform the psfgen operation to generate the PSF and PDB files for the bilayer patch from the packmol output.
        """
        self.next_basename(f'psfgen-{bilayer_name}')
        pg: PsfgenScripter = self.get_scripter('psfgen')
        pg.newscript(self.basename, additional_topologies=patch.addl_streamfiles)
        self.register_tops_streams_from_psfgen(pg.topologies)
        pg.usescript('bilayer_patch')
        pg.writescript(self.basename, guesscoord=False, regenerate=True, force_exit=True)
        state: StateArtifacts = self.get_current_artifact(f'{bilayer_name}_state')
        pdb: Path = state.pdb.path
        result = pg.runscript(pdb=pdb.name, o=self.basename)
        cell_to_xsc(patch.box, patch.origin, f'{self.basename}.xsc')
        patch.area = patch.box[0][0] * patch.box[1][1]
        self.register(dict(
                        psf=PSFFileArtifact(self.basename, pytestable=True), 
                        pdb=PDBFileArtifact(self.basename, pytestable=True), 
                        xsc=NAMDXscFileArtifact(self.basename)), 
                        key=f'{bilayer_name}_state', artifact_type=StateArtifacts)
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact)
        self.register_tmpfiles_from_tcl(f'{self.basename}.log')

    def pack_patch(self, patch: Bilayer, patch_name: str = None, seed=None, tolerance=None, nloop_all=200, nloop=200, half_mid_zgap=1.0, rotation_pm=20):
        """
        Packs the bilayer patch using Packmol.
        
        Parameters
        ----------
        seed : int, optional
            The random seed for the packing process. Default is None.
        tolerance : float, optional
            The tolerance for the packing process. Default is None.
        nloop_all : int, optional
            The total number of loops for the packing process. Default is 200.
        nloop : int, optional
            The number of loops for each individual structure in the packing process. Default is 200.
        half_mid_zgap : float, optional
            The half mid-plane gap in Å. Default is 1.0 Å.
        rotation_pm : float, optional
            The rotation angle in degrees for the patch. Default is 20.0 degrees.
        """
        self.next_basename(f'packmol-{patch_name}')
        pm: PackmolScripter = self.get_scripter('packmol')
        pm.newscript(self.basename)
        packmol_output_pdb = f'{self.basename}.pdb'
        pm.comment(f'packmol input automatically generated by pestifer {__pestifer_version__}')
        pm.addline(f'output {packmol_output_pdb}')
        pm.addline(f'filetype pdb')
        if seed is not None:
            pm.addline(f'seed {seed}')
        pm.addline(f'tolerance {tolerance}')
        pm.addline(f'nloop {nloop_all}')
        patch.write_packmol(pm, half_mid_zgap=half_mid_zgap, rotation_pm=rotation_pm, nloop=nloop)
        pm.writefile()
        result = pm.runscript()
        logger.debug(f'{self.basename} packmol result {result}')
        if result != 0:
            raise Exception(f'Packmol failed with result {result}')
        self.register(dict(pdb=PDBFileArtifact(self.basename, pytestable=True)), key=f'{patch_name}_state', artifact_type=StateArtifacts)
        if os.path.exists(f'{self.basename}.pdb_FORCED'):
            self.register(f'{self.basename}.pdb_FORCED', key=f'{patch_name}_packmol_forced', artifact_type=PackMolPDBForcedFileArtifact)
        self.register(self.basename, key='packmol_tcl', artifact_type=PackmolInputScriptArtifact)
        self.register(self.basename, key='packmol_log', artifact_type=PackmolLogFileArtifact)
        if os.path.exists(f'{self.basename}_packmol-results.yaml'):
            self.register(f'{self.basename}_packmol-results.yaml', key=f'{patch_name}_packmol_results', artifact_type=YAMLFileArtifact)
        csvs = glob.glob(f'{self.basename}*_packmol.csv')
        if len(csvs) > 0:
            logger.debug(f'CSVs:')
            my_logger(csvs, logger.debug)
            self.register([CSVDataFileArtifact(x) for x in csvs], key=f'{patch_name}_packmol_csvs', artifact_type=CSVDataFileArtifactList)
        pngs = glob.glob(f'{self.basename}*_packmol.png')
        if len(pngs) > 0:
            logger.debug(f'PNGs:')
            my_logger(pngs, logger.debug)
            self.register([PNGImageFileArtifact(x) for x in pngs], key=f'{patch_name}_packmol_pngs', artifact_type=PNGImageFileArtifactList)

    def equilibrate_bilayer(self, bilayer: Bilayer, bilayer_name: str, relaxation_protocol: list[dict] = None):
        """
        Equilibrates the bilayer patch using the specified user dictionary and relaxation protocol.
        
        Parameters
        ----------
        bilayer : Bilayer
            The bilayer object to be equilibrated.
        bilayer_name : str
            A string identifier of this bilayer used for keeping track of artifacts.
        relaxation_protocol : list, optional
            A list of dictionaries specifying the stages of the relaxation protocol.
            If not provided, a hard-coded relaxation protocol will be used. 
        """
        self.next_basename(f'equilibration-{bilayer_name}')
        # user_dict=deepcopy(self.config['user'])
        state: StateArtifacts = self.get_current_artifact(f'{bilayer_name}_state')
        logger.debug(f'Bilayer area before equilibration: {bilayer.area:.3f} {sA2_}')
        if not relaxation_protocol:
            logger.debug(f'Using hard-coded relaxation protocol for {self.basename}!!')
            relaxation_protocol=[
                {'md': dict(ensemble='minimize', minimize=1000)},
                {'md': dict(ensemble='NVT', nsteps=1000)},
                {'md': dict(ensemble='NPT', nsteps=200)},
                {'md': dict(ensemble='NPT', nsteps=400)},
                {'md': dict(ensemble='NPT', nsteps=800)},
                {'md': dict(ensemble='NPAT', nsteps=1600)},
                {'md': dict(ensemble='NPAT', nsteps=3200)},
                {'md': dict(ensemble='NPAT', nsteps=6400)},
                {'md': dict(ensemble='NPAT', nsteps=12800)},
                {'md': dict(ensemble='NPAT', nsteps=25600)}]
        else:
            logger.debug(f'Using user-specified relaxation protocol:')
            my_logger(relaxation_protocol, logger.debug)
        for stage in relaxation_protocol:
            specs = stage['md']
            specs['addl_paramfiles'] = bilayer.addl_streamfiles
            if specs.get('ensemble', None) in ['NPT', 'npt', 'NPAT', 'npat']:
                if not 'other_parameters' in specs:  # never true due to ycleptic base.yaml
                    specs['other_parameters'] = {'useflexiblecell': True, 'useconstantratio': True,
                                                 'pressureProfile': 'on', 'pressureProfileSlabs': 30,
                                                 'pressureProfileFreq': 100}
                else:
                    if not 'useflexiblecell' in specs['other_parameters']:
                        specs['other_parameters']['useflexiblecell'] = True
                    if not 'useconstantratio' in specs['other_parameters']:
                        specs['other_parameters']['useconstantratio'] = True
                    if self.provisions['processor-type'] != 'gpu':  # GPU NAMD 3.0.1 does not support pressure profiles
                        if not 'pressureProfile' in specs['other_parameters']:
                            specs['other_parameters']['pressureProfile'] = 'on'
                        if not 'pressureProfileSlabs' in specs['other_parameters']:
                            specs['other_parameters']['pressureProfileSlabs'] = 30
                        if not 'pressureProfileFreq' in specs['other_parameters']:
                            specs['other_parameters']['pressureProfileFreq'] = 100
        timeseries = ['density', ['a_x', 'b_y', 'c_z']]
        profiles = ['pressure']
        if self.provisions['processor-type'] != 'gpu':
            timeseries.append('pressure')  # To do: change this to pressureProfile plotting
        tasklist_user = [{'continuation': dict(psf=state.psf.name, pdb=state.pdb.name, xsc=state.xsc.name)}]
        tasklist_user.extend(relaxation_protocol)
        tasklist_user.extend([
            {'mdplot': dict(timeseries=timeseries, profiles=profiles, legend=True, grid=True, basename=self.basename)},
            # {'terminate': dict(basename=self.basename, cleanup=False, chainmapfile=f'{self.basename}-chainmap.yaml', statefile=f'{self.basename}-state.yaml')}
        ])
        subcontroller = self.subcontroller
        subcontroller.config['user']['title'] = f'Bilayer equilibration from {self.basename}'
        subcontroller.reconfigure_tasks(tasklist_user)
        for task in subcontroller.tasks:
            save_task_name = task.taskname
            task_name = f'{self.taskname}-{task.taskname}-{bilayer_name}'
            task.override_taskname(task_name)
            # logger.debug(f'Subcontroller overrides task name {save_task_name} with {task.taskname}')
        subcontroller.do_tasks()
        last_task = subcontroller.tasks[-1]
        bilayer_state: StateArtifacts = last_task.get_current_artifact('state')
        assert bilayer_state is not None
        assert bilayer_state.psf.exists()
        assert bilayer_state.pdb.exists()
        assert bilayer_state.vel.exists()
        assert bilayer_state.xsc.exists()
        assert bilayer_state.coor.exists()
        self.subcontroller.pipeline.rekey('state', f'{bilayer_name}_state')
        self.import_artifacts(subcontroller.pipeline)
        bilayer.box, bilayer.origin = cell_from_xsc(bilayer_state.xsc.name)
        bilayer.area = bilayer.box[0][0] * bilayer.box[1][1]
        logger.debug(f'{self.basename} area after equilibration: {bilayer.area:.3f} {sA2_}')

    def make_quilt_from_patch(self):
        """
        Create a quilt from the bilayer patch or patches.
        This method generates a quilt that combines the bilayer patches into a single structure.
        It uses the psfgen scripter to create a script that builds the quilt based on the provided patches.
        The quilt is then equilibrated and saved with the appropriate state variables.
        """
        
        logger.debug(f'Creating quilt from patch')
        self.next_basename('quilt')
        additional_topologies=[]
        if self.patch is not None:
            patch_state: StateArtifacts = self.get_current_artifact('patch_state')
            pdb = patch_state.pdb.name
            xsc = patch_state.xsc.name
            psf = patch_state.psf.name
            pdbA = pdbB = pdb
            psfA = psfB = psf
            xscA = xscB = xsc
        elif self.patchA is not None and self.patchB is not None:
            patchA_state: StateArtifacts = self.get_current_artifact('patchA_state')
            psfA = patchA_state.psf.name
            pdbA = patchA_state.pdb.name
            xscA = patchA_state.xsc.name
            patchB_state: StateArtifacts = self.get_current_artifact('patchB_state')
            psfB = patchB_state.psf.name
            pdbB = patchB_state.pdb.name
            xscB = patchB_state.xsc.name
        else:
            raise ValueError("No valid patch state found.")
        for patch in [self.patch, self.patchA, self.patchB]:
            if patch is None:
                continue
            additional_topologies += patch.addl_streamfiles
        additional_topologies = list(set(additional_topologies))
        pg: PsfgenScripter = self.get_scripter('psfgen')
        pg.newscript(self.basename, additional_topologies=additional_topologies)
        self.register_tops_streams_from_psfgen(pg.topologies)
        pg.usescript('bilayer_quilt')
        pg.writescript(self.basename, guesscoord=False, regenerate=False, force_exit=True, writepsf=False, writepdb=False)
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        margin = self.embed_specs.get('xydist', 10.0)
        protein_state: StateArtifacts = self.get_current_artifact('state')
        if protein_state is not None and protein_state.pdb.exists():
            # we will eventually embed a protein in here, so send its pdb along to help size the bilayer
            result = pg.runscript(propdb=protein_state.pdb.name, margin=margin, psfA=psfA, pdbA=pdbA, psfB=psfB, pdbB=pdbB, xscA=xscA, xscB=xscB, o=self.basename)
        else:
            dims = self.bilayer_specs.get('dims', [0,0])
            if dims is None or len(dims) != 2:
                dims = [0,0]
            dimx, dimy = dims
            npatch = self.bilayer_specs.get('npatch', [0,0])
            if npatch is None or len(npatch) != 2:
                npatch = [0,0]
            npatchx, npatchy = npatch
            if npatchx != 0 and npatchy != 0:
                result = pg.runscript(nx=npatchx, ny=npatchy, psfA=psfA, pdbA=pdbA,
                                      psfB=psfB, pdbB=pdbB, xscA=xscA, xscB=xscB, o=self.basename)
            elif dimx != 0 and dimy != 0:
                result = pg.runscript(dimx=dimx, dimy=dimy, psfA=psfA, pdbA=pdbA,
                                      psfB=psfB, pdbB=pdbB, xscA=xscA, xscB=xscB, o=self.basename)
        if result != 0:
            raise RuntimeError(f'psfgen failed with result {result} for {self.basename}')
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact)
        self.register_tmpfiles_from_tcl(f'{self.basename}.log')
        self.register(dict(
            pdb=PDBFileArtifact(self.basename, pytestable=True), 
            psf=PSFFileArtifact(self.basename, pytestable=True), 
            xsc=NAMDXscFileArtifact(self.basename)), key='quilt_state', artifact_type=StateArtifacts)
        self.quilt = Bilayer()
        self.quilt.addl_streamfiles = additional_topologies
        self.quilt.box, self.quilt.origin = cell_from_xsc(f'{self.basename}.xsc')
        self.quilt.area = self.quilt.box[0][0] * self.quilt.box[1][1]
        relaxation_protocol = self.bilayer_specs.get('relaxation_protocols', {}).get('quilt', {})
        self.equilibrate_bilayer(self.quilt,
                                 bilayer_name='quilt',
                                 relaxation_protocol=relaxation_protocol)

    def embed_protein(self):
        """
        Embed the protein into the bilayer patch or quilt.
        This method uses the psfgen scripter to create a script that embeds the protein
        into the bilayer based on the provided specifications.
        It retrieves the embedding specifications, such as head and tail groups, reference groups,
        and orientation settings, and runs the script to perform the embedding.
        The resulting state variables are updated with the new PSF, PDB, and coordinate files.
        """
        if not self.embed_specs:
            logger.debug('No embed specs.')
            return
        zvals = np.zeros(2)
        dum_pdb: PDBFileArtifact = self.get_current_artifact('base_coordinates_dum')
        if dum_pdb and dum_pdb.exists():
            logger.debug(f'Using DUM pdb {dum_pdb.name} for embedding coordinates')
            z_ref_group = ""
            dum_struct = PDBParser(filepath=dum_pdb.path).parse().parsed
            zvals = np.array(list(set(a.z for a in dum_struct['HETATM'])))
            logger.debug(f'DUM z-values: {zvals}')
            assert len(zvals) == 2, "DUM pdb must have exactly two distinct z-values for head and tail groups"
            z_value = zvals.mean()
        else:
            z_ref_group = self.embed_specs.get('z_ref_group', {}).get('text', None)
            z_value = self.embed_specs.get('z_ref_group', {}).get('z_value', 0.0)
        self.next_basename('embed')
        pg: PsfgenScripter = self.scripters['psfgen']
        pg.newscript(self.basename, additional_topologies=self.quilt.addl_streamfiles)
        pg.usescript('bilayer_embed')
        pg.writescript(self.basename, guesscoord=False, regenerate=True, force_exit=True, writepsf=False, writepdb=False)
        self.register(self.basename, key='tcl', artifact_type=PsfgenInputScriptArtifact)
        quilt_state: StateArtifacts = self.get_current_artifact('quilt_state')
        bilayer_psf: str = quilt_state.psf.name
        bilayer_pdb: str = quilt_state.pdb.name
        bilayer_xsc: str = quilt_state.xsc.name
        protein_state: StateArtifacts = self.get_current_artifact('state')
        protein_psf: str = protein_state.psf.name
        protein_pdb: str = protein_state.pdb.name
        logger.debug(f'Embedding {protein_pdb} with z_ref_group {z_ref_group} at z={z_value:.3f} into bilayer {bilayer_pdb}')
        result = pg.runscript(psf=protein_psf,
                              pdb=protein_pdb,
                              bilayer_psf=bilayer_psf,
                              bilayer_pdb=bilayer_pdb,
                              bilayer_xsc=bilayer_xsc,
                              z_ref_group=protect_str_arg(z_ref_group),
                              z_lo_dum=zvals[0],
                              z_hi_dum=zvals[1],
                              z_value=z_value,
                              o=self.basename)
        if result != 0:
            raise RuntimeError(f'psfgen failed with result {result} for {self.basename}')
        self.register(self.basename, key='log', artifact_type=PsfgenLogFileArtifact)
        self.register(dict(
            psf=PSFFileArtifact(self.basename, pytestable=True), 
            pdb=PDBFileArtifact(self.basename, pytestable=True), 
            xsc=NAMDXscFileArtifact(self.basename)), key='state', artifact_type=StateArtifacts)
        self.register_tmpfiles_from_tcl(f'{self.basename}.log')
        logger.debug(f'Embedding completed with result {result}')
        return result

    
