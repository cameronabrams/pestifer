# Author: Cameron F. Abrams
"""
Definition of the :class:`MDTask` class for handling molecular dynamics (MD) simulations using NAMD.
This class is a descendant of the :class:`BaseTask <pestifer.core.basetask.BaseTask>` class and is used to run NAMD simulations,
manage the state of the simulation, and handle various aspects of the MD task such as ensembles, constraints, and colvars.
It manages the setup and execution of NAMD runs, including reading and writing necessary files,
updating artifacts, and handling the results of the simulation.

Usage is described in the :ref:`subs_runtasks_md` documentation.
"""
import glob
import logging
import os

logger = logging.getLogger(__name__)

from .basetask import VMDTask
from ..core.artifacts import *

from ..scripters.namdscripter import NAMDScripter
from ..scripters.namdcolvarinputscripter import NAMDColvarInputScripter
from ..util.util import is_periodic

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
        logger.debug(f'namdrun called with baselabel {baselabel} and extras {extras}, script_only={script_only}')
        specs = self.specs
        addl_paramfiles: list[str] = specs.get('addl_paramfiles', [])
        logger.debug(f'md task specs {specs}')
        ensemble: str = specs['ensemble']
        if ensemble.casefold() not in ['NPAT'.casefold(), 'NPT'.casefold(), 'NVT'.casefold(), 'minimize'.casefold()]:
            raise Exception(f'Error: {ensemble} is not a valid ensemble type. Must be \'NPAT\', \'NPT\', \'NVT\', or \'minimize\'')
        if not baselabel:
            self.next_basename(ensemble)
        else:
            self.next_basename(baselabel)
        
        params = {}
        psf: Path = self.get_current_artifact_path('psf')
        pdb: Path = self.get_current_artifact_path('pdb')
        coor: Path = self.get_current_artifact_path('coor')
        vel: Path = self.get_current_artifact_path('vel')
        xsc: Path = self.get_current_artifact_path('xsc')
        prior_charmmff_parfiles: CharmmffParFileArtifacts = self.get_current_artifact('charmmff_parfiles')
        prior_paramfiles = []
        if prior_charmmff_parfiles:
            prior_paramfiles = [x.name for x in prior_charmmff_parfiles]
        firsttimestep = self.get_current_artifact_data('firsttimestep')
        if not firsttimestep:
            firsttimestep = 0
        if specs.get('firsttimestep', None) is not None:
            firsttimestep = specs['firsttimestep']
        if xsc is not None:
            self.register(DataArtifact(is_periodic(xsc.name), key='periodic'))

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
        params['structure'] = psf.name
        params['coordinates'] = pdb.name
        params['temperature'] = '$temperature'

        if coor:
            params['bincoordinates'] = coor.name
        if vel:
            params['binvelocities'] = vel.name
            del params['temperature']
        if xsc:
            params['extendedSystem'] = xsc.name
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
            writer.construct_on_pdb(specs=colvars, pdb=pdb)
            writer.writefile()
            self.register(NAMDColvarsConfigArtifact(f'{self.basename}-cv'))
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

        na: NAMDScripter = self.scripters['namd']
        na.newscript(self.basename, addl_paramfiles=list(set(addl_paramfiles + prior_paramfiles)))
        self.register(CharmmffParFileArtifacts([CharmmffParFileArtifact(x) for x in na.parameters if x.endswith('.prm')]), key='charmmff_parfiles')
        self.register(CharmmffStreamFileArtifacts([CharmmffStreamFileArtifact(x) for x in na.parameters if x.endswith('.str')]), key='charmmff_streamfiles')
        cpu_override = specs.get('cpu-override', False)
        logger.debug(f'CPU-override is {cpu_override}')
        na.writescript(params, cpu_override=cpu_override)
        self.register(NAMDConfigFileArtifact(self.basename))
        if not script_only:
            local_execution_only = not self.get_current_artifact_value('periodic')
            single_gpu_only = kwargs.get('single_gpu_only', False) or constraints
            result = na.runscript(single_molecule=local_execution_only, local_execution_only=local_execution_only, single_gpu_only=single_gpu_only, cpu_override=cpu_override)
        for at in [PDBFileArtifact, NAMDVelFileArtifact, NAMDXscFileArtifact, \
                   NAMDCoorFileArtifact, NAMDDcdFileArtifact, NAMDXstFileArtifact, \
                   NAMDLogFileArtifact, NAMDColvarsOutputArtifact]:
            artifact=at(self.basename)
            if artifact.exists():
                self.register(artifact)
            other_files = glob.glob(f'FFTW*txt')
            for f in other_files:
                if os.path.isfile(f):
                    self.register(TXTFileArtifact(f, key='NAMD FFTW warmup'))
            self.coor_to_pdb()
            if hasattr(na.logparser,'dataframes'):
                for key in na.logparser.dataframes:
                    artifact=CSVDataFileArtifact(f'{self.basename}-{key}')
                    if artifact.exists():
                        self.register(artifact, key=f'{key}-csv')
            if result != 0:
                return -1

            if ensemble.casefold() != 'minimize'.casefold():
                self.register(DataArtifact(firsttimestep+nsteps), key='firsttimestep')
            else:
                self.register(DataArtifact(firsttimestep+specs['minimize']), key='firsttimestep')

        return 0
