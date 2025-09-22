# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
A class for parsing the config file and creating the :class:`Config` object 
Pestifer's user-configuration input uses ycleptic, an enhanced, YAML-based
configuration file manager.  The :class:`Config` object is a descendent of the
:class:`Yclept` class.  It also houses the :class:`pestifer.core.resourcemanager.ResourceManager` object, which manages
access to the contents of Pestifer's :mod:`pestifer.resources` subpackage.
"""

import os
import logging
import shutil
from GPUtil import getGPUs
from importlib.metadata import version
from ycleptic.yclept import Yclept

from .resourcemanager import ResourceManager
from ..tasks.taskcollections import TaskList
from ..util.stringthings import my_logger
from ..scripters import PsfgenScripter, NAMDColvarInputScripter, PackmolScripter, VMDScripter, GenericScripter
from ..scripters.namdscripter import NAMDScripter

logger = logging.getLogger(__name__)

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
    def __init__(self, userfile='', userdict={}, quiet=False, RM: ResourceManager = None, basefile: str = ''):
        self.userfile = userfile
        self.userdict = userdict
        self.quiet = quiet
        self.basefile = basefile
        self.RM = RM
        self.my_processor_info = ''
        self.kwargs_to_scripters = {}

    def configure_new(self):
        if not self.quiet:
            vrep = f'pestifer v. {version("pestifer")}\n  ycleptic v. {version("ycleptic")}\n  pidibble v. {version("pidibble")}'
            my_logger(vrep, logger.info, just='<', frame='*', fill='', no_indent=True)
            # my_logger(str(self.RM), logger.debug, just='<', frame='*', fill='', no_indent=True)
        # resolve full pathname of YCleptic base config for this application

        self.RM = ResourceManager()
        return self.configure()
    
    def configure(self):
        if not self.basefile:
            self.basefile = self.RM.get_ycleptic_config()
        assert os.path.exists(self.basefile), f'Base config file {self.basefile} does not exist.'
        # ycleptic's init:
        super().__init__(self.basefile, userfile=self.userfile, userdict=self.userdict)

        self.my_processor_info = self._set_processor_info()
        if not self.quiet:
            my_logger(self.my_processor_info, logger.info, just='<', frame='*', fill='')
        self._set_internal_shortcuts()
        self._set_shell_commands(verify_access=(self.userfile != ''))
        self.RM.update_charmmff(
            tarball=self['user']['charmmff'].get('tarball', ''),
            user_custom_directory=self['user']['charmmff'].get('custom_directory', ''),
            user_pdbrepository_paths=self['user']['charmmff'].get('pdbrepository', [])
        )
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
        """ Create a taskless subconfiguration from the progenitor configuration. """
        subconfig = self.__class__(quiet=True, RM=self.RM).configure()
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
            namd_config = self['user']['namd'],
            namd_deprecates = self['user']['namd'].get('deprecated3', {}),
            namd_gpu = self.shell_commands['namd3gpu'],
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
        it defaults to the local CPU count. It also checks for available GPUs using GPUtil.
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
            gpus = getGPUs()
            if len(gpus) > 0:
                self.ngpus = len(gpus)
                self.gpu_devices = ','.join([str(x.id) for x in gpus])
                ess = 's' if self.ngpus > 1 else ''
                retstr += f'; {self.ngpus} gpu{ess}'

        self.ncpus = ncpus
        return retstr

    def _set_shell_commands(self, verify_access=True):
        """ Defines all shell commands used by Pestifer """
        required_commands = ['charmrun', 'namd3', 'vmd', 'catdcd', 'packmol']
        command_alternates = {'namd3': 'namd2'}
        self.shell_commands = {}
        for rq in required_commands:
            self.shell_commands[rq] = self['user']['paths'][rq]
            rq_resolved = shutil.which(self.shell_commands[rq])
            rq_alt = command_alternates.get(rq, None)
            if not rq_resolved and not rq_alt:
                raise FileNotFoundError(f'Cannot find required command {self.shell_commands[rq]} in your path.')
            
            if not rq_resolved:
                if rq in command_alternates:
                    rqalt = command_alternates[rq]
                    self.shell_commands[rq] = self['user']['paths'][rqalt]
                    altrq_resolved = shutil.which(self.shell_commands[rq])
                    if altrq_resolved:
                        logger.info(f'Using alternate command {self.shell_commands[rq]} for {rq}.')
                        rq_resolved = altrq_resolved
                    else:
                        raise FileNotFoundError(f'Cannot find required command {self.shell_commands[rq]} or alternate {self.shell_commands[rqalt]} in your path.')
                else:
                    raise FileNotFoundError(f'Cannot find required command {self.shell_commands[rq]} in your path.')
            if rq_resolved is not None and verify_access:
                assert os.access(rq_resolved, os.X_OK), f'You do not have permission to execute {rq_resolved}'
        self.namd_type = self['user']['namd']['processor-type']
        cn = 'namd3gpu'
        self.shell_commands[cn] = self['user']['paths'][cn]
        rq_resolved = shutil.which(self.shell_commands[cn])
        if not rq_resolved:
            logger.warning(f'{self.shell_commands[cn]}: not found')
        if self.shell_commands['namd3gpu'] == self.shell_commands['namd3']:
            logger.warning(f'CPU and GPU namd3 have the same path: {self.shell_commands["namd3gpu"]}')
        if self.namd_type == 'gpu': # make sure we can run the GPUResident namd3
            if verify_access:
                assert os.access(rq_resolved, os.X_OK), f'You do not have permission to execute {rq_resolved}'
        self.namd_deprecates = self['user']['namd']['deprecated3']

    def _set_internal_shortcuts(self):
        self.use_terminal_progress = len(self.slurmvars) == 0
        RM = self.RM
        self.tcl_root = RM.get_tcldir()
        assert os.path.exists(self.tcl_root)
        self.tcl_pkg_path = RM.get_tcl_pkgdir()
        assert os.path.exists(self.tcl_pkg_path)
        self.tcl_script_path = RM.get_tcl_scriptsdir()
        assert os.path.exists(self.tcl_script_path)
        self.vmd_startup_script = os.path.join(self.tcl_root, 'vmdrc.tcl')
        assert os.path.exists(self.vmd_startup_script)
        # self.charmmff_toppar_path=RM.get_charmmff_toppardir()
        # assert os.path.exists(self.charmmff_toppar_path)
        self.charmmff_custom_path = RM.get_charmmff_customdir()
        assert os.path.exists(self.charmmff_custom_path)
        self.user_charmmff_toppar_path = ''
        if hasattr(self, 'user'):
            self.user_charmmff_toppar_path = os.path.join(self['user']['charmmff'], 'toppar')
            assert os.path.exists(self.user_charmmff_toppar_path)
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
