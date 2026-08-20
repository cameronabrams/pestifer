# Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""
Definition of the :class:`MDTask` class for handling molecular dynamics (MD) simulations using NAMD.
This class is a descendant of the :class:`BaseTask <pestifer.tasks.basetask.BaseTask>` class and is used to run NAMD simulations,
manage the state of the simulation, and handle various aspects of the MD task such as ensembles, constraints, and colvars.
It manages the setup and execution of NAMD runs, including reading and writing necessary files,
updating artifacts, and handling the results of the simulation.

"""
import glob
import logging
import zlib


from .basetask import VMDTask
from ..core.artifacts import *
from ..core.errors import PestiferBuildError

from ..scripters.namd import NAMDScripter
from ..scripters.namd_colvar_input import NAMDColvarInputScripter
from ..util.density_convergence import total_atoms
from ..util.util import is_periodic

logger = logging.getLogger(__name__)

# 'NPgT' (preferred; 'npgt'/'NPGT' accepted) is the constant-lateral-ratio barostat
# ensemble pestifer historically called 'NPAT'.  Pestifer's 'NPAT' has never been true
# constant area -- it sets useFlexibleCell + useConstantRatio (x and y scale together,
# the cross-sectional area is free to change), not useConstantArea.  'NPAT' is still
# accepted but warns once that a future version will reserve it for genuine constant-area
# runs, at which point 'NPgT' will be the name for this constant-ratio behavior.
_VALID_ENSEMBLE_KEYS = ('minimize', 'nvt', 'npt', 'npgt')
_npat_deprecation_warned = False


def normalize_ensemble(ensemble: str) -> str:
    """Casefold an ensemble name and fold the legacy ``NPAT`` alias onto ``npgt``.

    Returns the canonical lowercase key (``minimize``/``nvt``/``npt``/``npgt``); raises for
    anything else.  Emits a one-time deprecation warning when ``NPAT`` is used."""
    global _npat_deprecation_warned
    ens = str(ensemble).casefold()
    if ens == 'npat':
        if not _npat_deprecation_warned:
            logger.warning(
                "ensemble 'NPAT' currently runs a constant-x:y-ratio barostat (useFlexibleCell "
                "+ useConstantRatio) -- this is NPgT, not constant area. Use 'NPgT' for this "
                "behavior; a future pestifer version will reserve 'NPAT' for true constant-area "
                "(useConstantArea) simulations.")
            _npat_deprecation_warned = True
        ens = 'npgt'
    if ens not in _VALID_ENSEMBLE_KEYS:
        raise PestiferBuildError(
            f"Error: {ensemble} is not a valid ensemble type. "
            f"Must be 'NPgT', 'NPAT', 'NPT', 'NVT', or 'minimize'")
    return ens

class MDTask(VMDTask):
    """ 
    A class for handling all NAMD runs
    """
    _yaml_header = 'md'
    """
    YAML header for the MDTask, used to identify the task in configuration files as part of a ``tasks`` list.
    """

    @classmethod
    def pipeline_contract(cls, specs):
        from .pipeline_contract import TaskContract, STATE, MD_OUTPUT
        # runs dynamics on the current system and produces timeseries for mdplot
        return TaskContract(requires=(STATE,), provides=(STATE, MD_OUTPUT))

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

    def record_md_run(self, ensemble: str, psf, csv_keys: list):
        """Append this NAMD run to the pipeline's chronological MD history.

        One entry per NAMD invocation, so a chunked task contributes one per chunk; entries
        sharing a ``(controller, task_index)`` are one logical stage.

        ``natoms`` is the point of the record.  A consumer needs to know which earlier dynamics
        describe *the system it is holding now*, and the answer is not "the tasks adjacent to me":
        a `solvate` between two MD stages makes the earlier one incommensurable, while a `validate`
        or a `ring_check` between them changes nothing.  Recording the atom count lets that be
        decided from the data instead of from task adjacency.
        """
        try:
            natoms = total_atoms(psf.path) if psf is not None else None
        except Exception as e:                      # a malformed/unreadable PSF must not kill the run
            logger.debug(f'{self.taskname}: could not count atoms for the MD history ({e})')
            natoms = None
        entry = dict(controller=getattr(self.pipeline, 'controller_index', 0),
                     task_index=self.index,
                     taskname=self.taskname,
                     basename=self.basename,
                     ensemble=ensemble,
                     natoms=natoms,
                     csv_keys=list(csv_keys),
                     csv_basename=self.basename)
        history = self.get_current_artifact('md_history')
        if history is None:
            self.register([entry], key='md_history', artifact_type=MDHistoryArtifact)
        else:
            # mutate in place and re-register: register() sees the same list object, so it stamps
            # the existing artifact rather than superseding it into the pipeline's own history
            history.data.append(entry)
            self.register(history.data, key='md_history', artifact_type=MDHistoryArtifact)
        logger.debug(f'{self.taskname}: MD history now holds {len(self.get_current_artifact("md_history").data)} run(s)')


    def _namd_seed(self):
        """This invocation's NAMD RNG seed, derived from the build-wide base seed.

        NAMD's own default seed comes from the clock, so two runs of one config produced
        different trajectories and neither could be reproduced.  Pestifer pins a base
        (``namd.seed``) and gives each NAMD invocation its own stream derived from it, so the
        build as a whole is reproducible while its individual runs stay independent of one
        another -- reusing one seed across every run in a build would correlate them.

        The base is hashed *together with* the stream key rather than XORed against it.  XOR
        leaves the base's low bits untouched, so a replica set seeded 27021972/3/4 -- which is
        exactly how the sweep scripts seed replicas -- produced NAMD seeds that differed only in
        their low bits (1112820144, 1112820145, 1112820146).  Adjacent seeds are not necessarily
        independent streams for a generator seeded by simple state initialization, and "the
        replicas were seeded one apart" is not a property worth having to defend.  Hashing the
        pair avalanches the base, so consecutive replicas get unrelated seeds.

        Returns 0 when the user has set ``namd.seed: 0``, meaning "leave it to NAMD".
        """
        base = int(self.namd_global_config.get('seed', 0) or 0)
        if not base:
            return 0
        stream = f'{base}:{self.controller_index}:{self.index}:{self.subtaskcount}'
        return (zlib.crc32(stream.encode()) % (2 ** 31 - 1)) or 1

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
        ens = normalize_ensemble(ensemble)   # casefold; 'NPAT' -> 'npgt' (warns once)
        if not baselabel:
            self.next_basename(ensemble)
        else:
            self.next_basename(baselabel)
        
        params = {}
        state: StateArtifacts = self.get_current_artifact('state')
        if not state.psf or not state.pdb:
            raise PestiferBuildError(f'No PSF or PDB file found for task {self.taskname}')
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
        if ens in ('npt', 'npgt'):
            pressure = specs['pressure']
        params['tcl'] = []
        params['tcl'].append(f'set temperature {temperature}')

        nsteps = specs['nsteps']
        dcdfreq = specs['dcdfreq']
        xstfreq = specs['xstfreq']

        constraints = specs.get('constraints',{})
        other_params = specs.get('other_parameters',{})
        colvars = specs.get('colvar_specs',{})
        ssrestraints = specs.get('ssrestraints',{})
        params.update(self.namd_global_config['generic'])
        seed = self._namd_seed()
        if seed:
            params['seed'] = seed

        # Report what this invocation runs, for run-record.json.  Accumulated rather than
        # assigned: a chunked equilibration calls namdrun repeatedly, and the record should show
        # the total the task ran, not the size of its last chunk.
        self.record_outcome(ensemble=ens, seed=seed or None,
                            temperature=float(temperature) if temperature is not None else None,
                            pressure=(float(pressure) if ens in ('npt', 'npgt') else None))
        self.outcome['steps'] = int(self.outcome.get('steps', 0)) + int(nsteps or 0)
        self.outcome['namd_runs'] = int(self.outcome.get('namd_runs', 0)) + 1
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
            if ens in ('npt', 'npgt') and xstfreq:
                params['xstfreq'] = xstfreq

        if self.get_current_artifact_data('periodic'):
            params.update(self.namd_global_config['solvated'])
        else:
            params.update(self.namd_global_config['vacuum'])

        if ens in ('npt', 'nvt', 'npgt'):
            params.update(self.namd_global_config['thermostat'])
            if ens == 'npt':
                if not self.get_current_artifact_data('periodic'):
                    raise PestiferBuildError(f'Cannot use barostat on a system without PBCs')
                params['tcl'].append(f'set pressure {pressure}')
                params.update(self.namd_global_config['barostat'])
            elif ens == 'npgt':
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
                raise PestiferBuildError(f'No constraint reference PDB file found for task {self.taskname}')
            params['constraints'] = 'on'
            params['consref'] = consref.name
            params['conskfile'] = consref.name
            params['conskcol'] = 'O'
        if colvars:
            writer: NAMDColvarInputScripter = self.scripters['namd_colvar']
            writer.newfile(f'{self.basename}-cv.in')
            writer.construct_on_pdb(specs=colvars, pdb=state.pdb.path)
            writer.writefile()
            self.register(f'{self.basename}-cv', key='in', artifact_type=NAMDColvarsConfigArtifact)
            params['colvars'] = 'on'
            params['colvarsconfig'] = self.get_current_artifact_path('in')
        if ssrestraints.get('enabled', False):
            eb_file = self.make_ssrestraints_file(ssrestraints)
            params['extraBonds'] = 'on'
            params['extraBondsFile'] = eb_file

        params.update(other_params)
        params.update(extras)
        params['firsttimestep'] = firsttimestep
        if ens == 'minimize':
            assert specs['minimize'] > 0, f'Error: you must specify how many minimization cycles'
            params['minimize'] = specs['minimize']
        elif ens in ('nvt', 'npt', 'npgt'):
            assert nsteps > 0, f'Error: you must specify how many time steps to run'
            params['run'] = nsteps

        skip_standard_params = kwargs.pop('skip_standard_params', False)
        na: NAMDScripter = self.get_scripter('namd')
        na.newscript(self.basename, addl_paramfiles=list(set(addl_paramfiles + paramfilenames)),
                     skip_standard_params=skip_standard_params)
        # update the list of parfiles and streamfiles in the pipeline
        charmmff_parfiles.extend(x for x in na.parameters if x.endswith('.prm') and x not in paramfilenames)
        charmmff_streamfiles.extend(x for x in na.parameters if x.endswith('.str') and x not in paramfilenames)

        cpu_override = specs.get('cpu-override', False)
        logger.debug(f'CPU-override is {cpu_override}')
        na.writescript(params, cpu_override=cpu_override)
        logger.debug(f'Consolidating CHARMMFF parameters for {self.basename} using PSF {state.psf.name}')
        minimal_prm = na.consolidate_params(state.psf.name)
        if minimal_prm:
            logger.debug(f'Parameter consolidation complete; registering {minimal_prm} as artifact')
            self.register(minimal_prm, key='charmmff_minimal_prm', artifact_type=CharmmffParFileArtifact)
        else:
            logger.warning(f'Parameter consolidation skipped for {self.basename}; using full parameter set')
        self.register(self.basename, key='namd', artifact_type=NAMDConfigFileArtifact, pytestable=True)
        result = 0  # anticipate success
        if script_only:
            logger.debug(f'namdrun script_only is True; skipping execution')
            # self.pipeline.show_artifacts(header=f'Current artifacts after preparing NAMD run script only for {self.basename}')
            return result
        local_execution_only = not self.get_current_artifact_data('periodic')
        single_gpu_only = kwargs.get('single_gpu_only', False) or constraints
        single_cpu_only = specs.get('single-core', False)
        result = na.runscript(single_molecule=local_execution_only, local_execution_only=local_execution_only, single_gpu_only=single_gpu_only, cpu_override=cpu_override, single_cpu_only=single_cpu_only)
        if result != 0:
            raise PestiferBuildError(f'md task {self.taskname} failed.')
        coor = NAMDCoorFileArtifact(self.basename)
        vel = NAMDVelFileArtifact(self.basename)
        dcd = NAMDDcdFileArtifact(self.basename)
        if not coor.exists() and not vel.exists() and not dcd.exists():
            raise PestiferBuildError(f'NAMD run {self.basename} failed: no .coor, .vel, or .dcd files produced')
        xsc = NAMDXscFileArtifact(self.basename)
        xst = NAMDXstFileArtifact(self.basename)
        log = NAMDLogFileArtifact(self.basename)
        cvtraj = NAMDColvarsTrajectoryArtifact(self.basename)
        cvstate = NAMDColvarsStateArtifact(self.basename)
        
        # Restore chain IDs from the pre-MD state PDB: catdcd re-derives them from segid,
        # which would collapse any chain assignment that differs from the segid (e.g. merged
        # copies given distinct chains). state.pdb is guaranteed present (asserted above).
        self.coor_to_pdb(f'{self.basename}.coor', state.psf.name, reference_pdb=state.pdb.name)
        pdb = PDBFileArtifact(f'{self.basename}.pdb', pytestable=True)
        psf = state.psf
        if not pdb.exists():
            raise PestiferBuildError(f'NAMD run {self.basename} failed: no .pdb file produced from .coor')
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
        csv_keys = []
        if hasattr(na.logparser, 'dataframes'):
            for key in na.logparser.dataframes:
                self.register(f'{self.basename}-{key}', artifact_type=CSVDataFileArtifact, key=f'{key}-csv')
                csv_keys.append(key)
        self.record_md_run(ensemble=ens, psf=state.psf, csv_keys=csv_keys)
        if result != 0:
            return -1

        if ens != 'minimize':
            self.register(firsttimestep+nsteps, key='firsttimestep')
        else:
            self.register(firsttimestep+specs['minimize'], key='firsttimestep')

        self.pipeline.show_artifacts(header=f'Current artifacts after completing NAMD run {self.basename}')
        return 0
