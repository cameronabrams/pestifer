# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`MDTask` class for handling molecular dynamics (MD) simulations using NAMD.
This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and is used to run NAMD simulations,
manage the state of the simulation, and handle various aspects of the MD task such as ensembles, constraints, and colvars.
It manages the setup and execution of NAMD runs, including reading and writing necessary files,
updating artifacts, and handling the results of the simulation.

"""
import glob
import logging


from .basetask import VMDTask
from ..core.artifacts import *

from ..scripters.namdscripter import NAMDScripter
from ..scripters.namdcolvarinputscripter import NAMDColvarInputScripter
from ..util.util import is_periodic

logger = logging.getLogger(__name__)

class MDTask(VMDTask):
    """ 
    A class for handling all NAMD runs
    """
    _yaml_header = 'md'
    """
    YAML header for the MDTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """
    
    def provision(self, packet: dict):
        """
        Provision the MDTask with the provided packet of data.
        This method updates the task's provisions with the necessary configuration for NAMD runs.
        It sets up the global NAMD configuration and prepares the task for execution.
        
        Parameters
        ----------
        packet : dict
            A dictionary containing the necessary configuration data for the MDTask.
        """
        super().provision(packet)
        ensemble = self.specs.get('ensemble', None)
        self.extra_message = ''
        if ensemble:
            self.extra_message = f'ensemble: {ensemble}'
        self.namd_global_config: dict = self.provisions.get('namd_global_config', None)

    def do(self):
        """
        Execute the MD task; special messaging is handled here, and do() is stubbed out.
        """
        self.result = self.namdrun()
        return self.result

    def namdrun(self, baselabel='', extras={}, script_only=False, **kwargs):
        """
        Run a NAMD simulation based on the specifications provided in the task.
        This method prepares the necessary parameters, writes the NAMD script, and executes it.
        It handles different ensembles (NPT, NVT, minimize), sets up constraints and colvars, and manages the simulation state.
        If `script_only` is True, it only writes the script without executing it.
        
        Parameters
        ----------
        baselabel : str, optional
            The base label for the NAMD run. If not provided, it will be generated based on the ensemble type.
        extras : dict, optional
            Additional parameters to be included in the NAMD script.
        script_only : bool, optional
            If True, only the script will be written without executing it.
        **kwargs : dict, optional
            Additional keyword arguments that may include execution options such as `single_gpu_only`.
        """
        # logger.debug(f'namdrun called with baselabel {baselabel} and extras {extras}, script_only={script_only}')
        specs = self.specs
        addl_paramfiles: list[str] = specs.get('addl_paramfiles', [])
        logger.debug(f'MD task specs:')
        my_logger(specs, logger.debug)
        ensemble: str = specs['ensemble']
        if ensemble.casefold() not in ['NPAT'.casefold(), 'NPT'.casefold(), 'NVT'.casefold(), 'minimize'.casefold()]:
            raise Exception(f'Error: {ensemble} is not a valid ensemble type. Must be \'NPAT\', \'NPT\', \'NVT\', or \'minimize\'')
        if not baselabel:
            self.next_basename(ensemble)
        else:
            self.next_basename(baselabel)
        
        params = {}
        state: StateArtifacts = self.get_current_artifact('state')
        if not state.psf or not state.pdb:
            raise RuntimeError(f'No PSF or PDB file found for task {self.taskname}')
        charmmff_parfiles: CharmmffParFileArtifacts = self.get_current_artifact('charmmff_parfiles')
        if charmmff_parfiles is None:
            charmmff_parfiles = self.register([], key='charmmff_parfiles', artifact_type=CharmmffParFileArtifacts)
        charmmff_streamfiles: CharmmffStreamFileArtifacts = self.get_current_artifact('charmmff_streamfiles')
        if charmmff_streamfiles is None:
            charmmff_streamfiles = self.register([], key='charmmff_streamfiles', artifact_type=CharmmffStreamFileArtifacts)
        paramfilenames = [x.name for x in charmmff_parfiles] if charmmff_parfiles else []
        paramfilenames.extend([x.name for x in charmmff_streamfiles] if charmmff_streamfiles else [])
        firsttimestep = self.get_current_artifact_data('firsttimestep')
        if not firsttimestep:
            firsttimestep = 0
        if specs.get('firsttimestep', None) is not None:
            firsttimestep = specs['firsttimestep']
        if state.xsc is not None:
            self.register(is_periodic(state.xsc.name), key='periodic', artifact_type=DataArtifact)
        xsc = state.xsc
        coor = state.coor
        temperature = specs['temperature']
        if ensemble.casefold() == 'NPT'.casefold() or ensemble.casefold() == 'NPAT'.casefold():
            pressure = specs['pressure']
        params['tcl'] = []
        params['tcl'].append(f'set temperature {temperature}')

        nsteps = specs['nsteps']
        dcdfreq = specs['dcdfreq']
        xstfreq = specs['xstfreq']

        constraints = specs.get('constraints',{})
        other_params = specs.get('other_parameters',{})
        colvars = specs.get('colvar_specs',{})
        params.update(self.namd_global_config['generic'])
        params['structure'] = state.psf.name
        params['coordinates'] = state.pdb.name
        params['temperature'] = '$temperature'

        if state.coor:
            params['bincoordinates'] = state.coor.name
        if state.vel:
            params['binvelocities'] = state.vel.name
            del params['temperature']
        if state.xsc:
            params['extendedSystem'] = state.xsc.name
            if (ensemble.casefold()=='NPT'.casefold() or ensemble.casefold()=='NPAT'.casefold()) and xstfreq:
                params['xstfreq'] = xstfreq

        if self.get_current_artifact_data('periodic'):
            params.update(self.namd_global_config['solvated'])
        else:
            params.update(self.namd_global_config['vacuum'])

        if ensemble.casefold() in ['NPT'.casefold(),'NVT'.casefold(), 'NPAT'.casefold()]:
            params.update(self.namd_global_config['thermostat'])
            if ensemble.casefold() == 'NPT'.casefold():
                if not self.get_current_artifact_data('periodic'):
                    raise Exception(f'Cannot use barostat on a system without PBCs')
                params['tcl'].append(f'set pressure {pressure}')
                params.update(self.namd_global_config['barostat'])
            elif ensemble.casefold() == 'NPAT'.casefold():
                params['tcl'].append(f'set pressure {pressure}')
                params.update(self.namd_global_config['membrane'])
        params['outputName'] = f'{self.basename}'

        if dcdfreq:
            params['dcdfreq'] = dcdfreq
        if 'restartfreq' in specs:
            params['restartfreq'] = specs['restartfreq']
        if constraints:
            self.make_constraint_pdb(constraints)
            consref = self.get_current_artifact_path('consref')
            if not consref:
                raise RuntimeError(f'No constraint reference PDB file found for task {self.taskname}')
            params['constraints'] = 'on'
            params['consref'] = consref.name
            params['conskfile'] = consref.name
            params['conskcol'] = 'O'
        if colvars:
            writer: NAMDColvarInputScripter = self.scripters['namd_colvar']
            writer.newfile(f'{self.basename}-cv.in')
            writer.construct_on_pdb(specs=colvars, pdb=state.pdb)
            writer.writefile()
            self.register('{self.basename}-cv', key='in', artifact_type=NAMDColvarsConfigArtifact)
            params['colvars'] = 'on'
            params['colvarsconfig'] = self.get_current_artifact_path('in')

        params.update(other_params)
        params.update(extras)
        params['firsttimestep'] = firsttimestep
        if ensemble.casefold() == 'minimize'.casefold():
            assert specs['minimize'] > 0, f'Error: you must specify how many minimization cycles'
            params['minimize'] = specs['minimize']
        elif ensemble.casefold() in ['NVT'.casefold(), 'NPT'.casefold(), 'NPAT'.casefold()]:
            assert nsteps > 0, f'Error: you must specify how many time steps to run'
            params['run'] = nsteps

        na: NAMDScripter = self.get_scripter('namd')
        na.newscript(self.basename, addl_paramfiles=list(set(addl_paramfiles + paramfilenames)))
        # update the list of parfiles and streamfiles in the pipeline
        charmmff_parfiles.extend(x for x in na.parameters if x.endswith('.prm') and x not in paramfilenames)
        charmmff_streamfiles.extend(x for x in na.parameters if x.endswith('.str') and x not in paramfilenames)

        cpu_override = specs.get('cpu-override', False)
        logger.debug(f'CPU-override is {cpu_override}')
        na.writescript(params, cpu_override=cpu_override)
        self.register(self.basename, key='namd', artifact_type=NAMDConfigFileArtifact, pytestable=True)
        result = 0  # anticipate success
        if script_only:
            logger.debug(f'namdrun script_only is True; skipping execution')
            # self.pipeline.show_artifacts(header=f'Current artifacts after preparing NAMD run script only for {self.basename}')
            return result
        local_execution_only = not self.get_current_artifact_data('periodic')
        single_gpu_only = kwargs.get('single_gpu_only', False) or constraints
        result = na.runscript(single_molecule=local_execution_only, local_execution_only=local_execution_only, single_gpu_only=single_gpu_only, cpu_override=cpu_override)
        if result != 0:
            raise RuntimeError(f'md task {self.taskname} failed.')
        coor = NAMDCoorFileArtifact(self.basename)
        vel = NAMDVelFileArtifact(self.basename)
        dcd = NAMDDcdFileArtifact(self.basename)
        if not coor.exists() and not vel.exists() and not dcd.exists():
            raise RuntimeError(f'NAMD run {self.basename} failed: no .coor, .vel, or .dcd files produced')
        xsc = NAMDXscFileArtifact(self.basename)
        xst = NAMDXstFileArtifact(self.basename)
        log = NAMDLogFileArtifact(self.basename)
        cvtraj = NAMDColvarsTrajectoryArtifact(self.basename)
        cvstate = NAMDColvarsStateArtifact(self.basename)
        
        self.coor_to_pdb(f'{self.basename}.coor', state.psf.name)
        pdb = PDBFileArtifact(f'{self.basename}.pdb', pytestable=True)
        psf = state.psf
        if not pdb.exists():
            raise RuntimeError(f'NAMD run {self.basename} failed: no .pdb file produced from .coor')
        self.register(dict(pdb=pdb,
                           coor=coor,
                           vel=vel,
                           xsc=xsc if xsc.exists() else None,
                           psf=psf), key='state', artifact_type=StateArtifacts)  # an md run cannot change the PSF file
        falist = [x for x in [dcd, xst, log, cvtraj, cvstate] if x and x.exists()]
        self.register(falist, artifact_type=FileArtifactList, key='md_output_files')
        other_files = glob.glob(f'FFTW*txt')
        for f in other_files:
            self.register(f, artifact_type=TXTFileArtifact, key='NAMD FFTW warmup')
        if hasattr(na.logparser, 'dataframes'):
            for key in na.logparser.dataframes:
                self.register(f'{self.basename}-{key}', artifact_type=CSVDataFileArtifact, key=f'{key}-csv')
        if result != 0:
            return -1

        if ensemble.casefold() != 'minimize'.casefold():
            self.register(firsttimestep+nsteps, key='firsttimestep')
        else:
            self.register(firsttimestep+specs['minimize'], key='firsttimestep')

        self.pipeline.show_artifacts(header=f'Current artifacts after completing NAMD run {self.basename}')
        return 0
