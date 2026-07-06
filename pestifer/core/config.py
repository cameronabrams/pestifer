# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
A class for parsing the config file and creating the :class:`Config` object 
Pestifer's user-configuration input uses ycleptic, an enhanced, YAML-based
configuration file manager.  The :class:`Config` object is a descendent of the
:class:`Yclept` class.  It also houses the :class:`pestifer.core.resourcemanager.ResourceManager` object, which manages
access to the contents of Pestifer's :mod:`pestifer.resources` subpackage.
"""

import os
import sys
import copy
import logging
import shutil
import subprocess
from importlib.metadata import version
from importlib.resources import files as pkg_files
from ycleptic import Yclept

from .errors import PestiferError
from .resourcemanager import ResourceManager
from ..tasks.taskcollections import TaskList
from ..util.stringthings import my_logger
from ..scripters import PsfgenScripter, NAMDColvarInputScripter, PackmolScripter, VMDScripter, GenericScripter
from ..scripters.namd import NAMDScripter
from ..scripters.packmol import check_packmol_version

logger = logging.getLogger(__name__)


def detect_local_gpu_ids() -> list:
    """Return the indices of local NVIDIA GPUs by querying ``nvidia-smi``.

    Returns an empty list if ``nvidia-smi`` is not on PATH, fails, or reports no
    GPUs.  This replaces the former GPUtil dependency, which was an unmaintained
    wrapper around the same ``nvidia-smi`` query and pulled in ``distutils`` (and
    hence ``setuptools``) on Python 3.12+.
    """
    exe = shutil.which('nvidia-smi')
    if not exe:
        return []
    try:
        out = subprocess.run([exe, '--query-gpu=index', '--format=csv,noheader'],
                             capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return []
    if out.returncode != 0:
        return []
    return [int(tok) for tok in out.stdout.split() if tok.strip().isdigit()]


class Config(Yclept):
    """ 
    A class for managing the configuration of Pestifer.
    This class extends the :class:`ycleptic.yclept.Yclept` class to provide additional functionality
    specific to Pestifer's configuration needs.

    Parameters
    ----------
    userfile : str, optional
        Path to the user-specific configuration file. If not provided, the default configuration is used.
    userdict : dict, optional
        A dictionary of user-specific configuration options. If not provided, the default configuration is used.
    quiet : bool, optional
        If True, suppresses output to the console. Default is False.
    RM : ResourceManager
        An instance of the ResourceManager class, which manages access to the contents of Pestifer's resources.
    basefile : str
        Optional name of the Ycleptic-format base file
    """
    def __init__(self, userfile='', userdict={}, quiet=False, RM: ResourceManager = None, basefile: str = '', ncpus_override: int = 0):
        self.userfile = userfile
        self.userdict = userdict
        self.quiet = quiet
        self.basefile = basefile
        self.RM = RM
        self.ncpus_override = ncpus_override
        self.my_processor_info = ''
        self.kwargs_to_scripters = {}

    def configure_new(self):
        if not self.quiet:
            vrep = f'pestifer v. {version("pestifer")}\n  ycleptic v. {version("ycleptic")}\n  pidibble v. {version("pidibble")}'
            my_logger(vrep, logger.info, just='<', frame='*', fill='', no_indent=True)
        return self.configure()

    def configure(self):
        if not self.basefile:
            self.basefile = str(pkg_files('pestifer.schema').joinpath('base.yaml'))
        assert os.path.exists(self.basefile), f'Base config file {self.basefile} does not exist.'
        # ycleptic's init:
        super().__init__(self.basefile, userfile=self.userfile, userdict=self.userdict)

        self.my_processor_info = self._set_processor_info()
        if not self.quiet:
            my_logger(self.my_processor_info, logger.info, just='<', frame='*', fill='')
        if self.RM is None:
            self.RM = ResourceManager(charmmff_config=self['user']['charmmff'])
        self._set_internal_shortcuts()
        self._set_shell_commands(verify_access=(self.userfile != ''))
        self._set_kwargs_to_scripters()
        self.scripters = {
            'psfgen': PsfgenScripter(**self.kwargs_to_scripters),
            'namd': NAMDScripter(**self.kwargs_to_scripters),
            'packmol': PackmolScripter(**self.kwargs_to_scripters),
            'tcl': VMDScripter(**self.kwargs_to_scripters),
            'data': GenericScripter(**self.kwargs_to_scripters),
            'vmd': VMDScripter(**self.kwargs_to_scripters),
            'namd_colvar': NAMDColvarInputScripter(**self.kwargs_to_scripters)
        }
        return self

    def taskless_subconfig(self) -> 'Config':
        """ Create a taskless subconfiguration from the progenitor configuration.

        The subconfiguration inherits the progenitor's user settings (NAMD launcher,
        paths, force field, psfgen options, ...) so that tasks executed by a
        subcontroller -- e.g. the relaxation MD inside ``make_membrane_system`` -- run
        in the same execution environment as the top-level run.  Only the task list is
        reset to empty; without this inheritance the subcontroller would fall back to
        schema defaults (e.g. ``cpu-parallel-launcher: auto``) and ignore user overrides.
        """
        user_overrides = {k: copy.deepcopy(v) for k, v in self['user'].items() if k != 'tasks'}
        subconfig = self.__class__(userdict=user_overrides, quiet=True, RM=self.RM).configure()
        subconfig['user']['tasks'] = TaskList([])
        return subconfig

    def get_scripter(self, scripter_name: str):
        """ 
        Get a scripter instance by name.

        Parameters
        ----------
        scripter_name : str
            The name of the scripter to retrieve. Must be one of 'psfgen', 'namd', 'packmol', 'tcl', 'data', or 'vmd'.

        Returns
        -------
        Scripter
            An instance of the requested scripter.
        """
        if scripter_name in self.scripters:
            return self.scripters[scripter_name]
        else:
            raise ValueError(f'Scripter {scripter_name} not found.')

    def _set_kwargs_to_scripters(self):
        """
        Defines a collection of resources sent to all scripters
        """
        assert self.RM.charmmff_content is not None, 'ResourceManager must have charmmff_content set.'
        self.kwargs_to_scripters = dict(
            charmmff = self.RM.charmmff_content,
            charmmff_content = self.RM.charmmff_content,
            charmmff_config = self['user']['charmmff'],
            charmrun = self.shell_commands['charmrun'],
            gpu_devices = self.gpu_devices,
            local_ncpus = self.local_ncpus,
            namd = self.shell_commands['namd3'],
            namd3gpu = self.shell_commands['namd3gpu'],
            namd_config = self['user']['namd'],
            namd_deprecates = self['user']['namd'].get('deprecated3', {}),
            namd_type = self.namd_type,
            namd_version = self['user']['namd'].get('version', '3.0'),
            ncpus = self.ncpus,
            ngpus = self.ngpus,
            packmol = self.shell_commands['packmol'],
            psfgen_config = self['user']['psfgen'],
            progress = self.use_terminal_progress,
            slurmvars = self.slurmvars,
            tcl_root = self.tcl_root,
            tcl_pkg_path = self.tcl_pkg_path,
            tcl_script_path = self.tcl_script_path,
            vmd = self.shell_commands['vmd'],
            vmd_startup_script = self.vmd_startup_script,
        )

    def _set_processor_info(self):
        """ 
        Determine the number of CPUs and GPUs available for this process.
        This method checks if the process is running under SLURM (a job scheduler)
        and retrieves the number of nodes and tasks per node. If SLURM variables are not set,
        it defaults to the local CPU count. It also checks for available GPUs via nvidia-smi.
        If running under SLURM, the return value includes the number of nodes and tasks per node.
        If running locally, it includes the local CPU count and the number of GPUs detected.
        If no GPUs are detected, it will not mention GPUs in the output.
        If GPUs are detected, it will include the number of GPUs and their IDs.

        Returns
        -------
        str
            A string summarizing the number of CPUs and GPUs available.
        """
        self.slurmvars = {k: os.environ[k] for k in os.environ if 'SLURM' in k}
        self.local_ncpus = os.cpu_count()
        self.gpus_allocated = ''
        self.ngpus = 0
        self.gpu_devices = ''
        retstr = ''
        if self.slurmvars and 'SLURM_NNODES' in self.slurmvars and 'SLURM_NTASKS_PER_NODE' in self.slurmvars:
            # we are in a batch execution managed by slurm
            nnodes = int(self.slurmvars['SLURM_NNODES'])
            ntaskspernode = int(self.slurmvars['SLURM_NTASKS_PER_NODE'])
            ncpus = nnodes * ntaskspernode
            retstr += f'SLURM: {nnodes} nodes; {ncpus} cpus'
            if 'SLURM_JOB_GPUS' in self.slurmvars:
                self.gpu_devices = self.slurmvars['SLURM_JOB_GPUS']
                self.ngpus = len(self.gpus_allocated.split(','))
                ess = 's' if self.ngpus > 1 else ''
                retstr += f'; {self.ngpus} gpu{ess}'
        else:
            retstr += f'Local: {self.local_ncpus} cpus'
            ncpus = self.local_ncpus
            gpu_ids = detect_local_gpu_ids()
            if gpu_ids:
                self.ngpus = len(gpu_ids)
                self.gpu_devices = ','.join(str(i) for i in gpu_ids)
                ess = 's' if self.ngpus > 1 else ''
                retstr += f'; {self.ngpus} gpu{ess}'

        if self.ncpus_override > 0:
            ncpus = self.ncpus_override
            retstr += f'; will use {ncpus} PEs (--ncpus override)'
        else:
            user_ncpus = self['user']['namd'].get('ncpus', 0)
            if user_ncpus > 0:
                ncpus = user_ncpus
                retstr += f'; will use {ncpus} PEs (config-specified)'
            else:
                retstr += f'; will use {ncpus} PEs (auto-detected)'
        self.ncpus = ncpus
        return retstr

    def _set_shell_commands(self, verify_access=True):
        """ Defines all shell commands used by Pestifer.

        When ``verify_access`` is True -- a build run, i.e. a user config file was
        supplied -- every command in ``required_commands`` must resolve on PATH and
        be executable, or a :class:`PestiferError` is raised.  When False -- a
        standalone post-processing subcommand such as ``desolvate`` or ``mdplot``,
        which uses only a subset (typically ``vmd``/``catdcd``) -- a missing command
        is recorded and warned about rather than treated as fatal, so those
        subcommands are not gated on the full NAMD/packmol toolchain being loaded.
        The task that actually needs a given command fails loudly itself if the
        command is missing when used (see e.g. ``DesolvateTask``).
        """
        required_commands = ['charmrun', 'namd3', 'vmd', 'catdcd', 'packmol']
        command_alternates = {'namd3': 'namd2'}
        self.shell_commands = {}

        def _missing(msg):
            if verify_access:
                raise PestiferError(msg)
            logger.warning(f'{msg} Continuing; the task that uses it will fail if it is actually needed.')

        for rq in required_commands:
            self.shell_commands[rq] = self['user']['paths'][rq]
            rq_resolved = shutil.which(self.shell_commands[rq])
            if not rq_resolved and rq in command_alternates:
                rqalt = command_alternates[rq]
                self.shell_commands[rq] = self['user']['paths'][rqalt]
                rq_resolved = shutil.which(self.shell_commands[rq])
                if rq_resolved:
                    logger.info(f'Using alternate command {self.shell_commands[rq]} for {rq}.')
            if not rq_resolved:
                if rq in command_alternates:
                    _missing(f'Cannot find or execute required command '
                             f'{self["user"]["paths"][rq]!r} or alternate '
                             f'{self["user"]["paths"][command_alternates[rq]]!r}.')
                else:
                    _missing(f'Cannot find or execute required command {self.shell_commands[rq]!r}.')
            if rq_resolved is not None and verify_access:
                assert os.access(rq_resolved, os.X_OK), f'You do not have permission to execute {rq_resolved}'
        if verify_access:
            try:
                check_packmol_version(self.shell_commands['packmol'])
            except RuntimeError as e:
                raise PestiferError(str(e)) from e
        namd3_path = self.shell_commands['namd3']
        namd3gpu_path = self['user']['paths']['namd3gpu']
        self.shell_commands['namd3gpu'] = namd3gpu_path
        if namd3gpu_path != namd3_path:
            namd3gpu_resolved = shutil.which(namd3gpu_path)
            if namd3gpu_resolved:
                if verify_access:
                    assert os.access(namd3gpu_resolved, os.X_OK), f'You do not have permission to execute {namd3gpu_resolved}'
                self.namd_type = 'gpu'
            else:
                logger.warning(f'namd3gpu {namd3gpu_path!r} not found; falling back to CPU mode')
                self.shell_commands['namd3gpu'] = namd3_path
                self.namd_type = 'cpu'
        else:
            self.namd_type = 'cpu'
        self.namd_deprecates = self['user']['namd']['deprecated3']

    def _set_internal_shortcuts(self):
        # Progress bars are a terminal animation (carriage-return redraws) rendered by
        # progressbar2 on stderr -- the same stream pestifer logs to. When stderr is not a
        # TTY (redirected to a log file, e.g. 'pestifer run x.yaml &> run.log'), each redraw
        # becomes a separate line and floods the log, so only enable them for an interactive
        # terminal session (and never under a batch scheduler).
        self.use_terminal_progress = len(self.slurmvars) == 0 and sys.stderr.isatty()
        RM = self.RM
        self.tcl_root = RM.get_tcldir()
        assert os.path.exists(self.tcl_root)
        self.tcl_pkg_path = RM.get_tcl_pkgdir()
        assert os.path.exists(self.tcl_pkg_path)
        self.tcl_script_path = RM.get_tcl_scriptsdir()
        assert os.path.exists(self.tcl_script_path)
        self.vmd_startup_script = os.path.join(self.tcl_root, 'vmdrc.tcl')
        assert os.path.exists(self.vmd_startup_script)
        self.namd_config_defaults = self['user']['namd']

        self.segtypes = RM.labels.segtypes
        if self['user']['psfgen']['segtypes']:
            RM.labels.update_segtypes(self['user']['psfgen']['segtypes'])
        # logger.debug(f'segtypes: {self.segtypes}')
        for atom_alias in RM.labels.aliases['atom']:
            if atom_alias not in self['user']['psfgen']['aliases']['atom']:
                # add the atom alias to the user config
                self['user']['psfgen']['aliases']['atom'].append(atom_alias)
        for residue_alias in RM.labels.aliases['residue']:
            if residue_alias not in self['user']['psfgen']['aliases']['residue']:
                # add the residue alias to the user config
                self['user']['psfgen']['aliases']['residue'].append(residue_alias)
        RM.labels.update_aliases(residue_aliases=self['user']['psfgen']['aliases']['residue'],
                                 atom_aliases=self['user']['psfgen']['aliases']['atom'])
        # logger.debug(f'psfgen aliases: {self["user"]["psfgen"]["aliases"]}')
