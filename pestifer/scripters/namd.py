# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import logging
import os

from .tclscripter import TcLScripter
from ..core.command import Command
from ..logparsers import NAMDLogParser, NAMDxstParser
from ..util.progress import NAMDProgress

logger = logging.getLogger(__name__)

class NAMDScripter(TcLScripter):
    """
    This class extends the TcLScripter class to provide functionality for creating and managing NAMD scripts (which NAMD refers to as "configurations" -- this is not to be confused with Pestifer's Config class).

    Parameters
    ----------
    config : Config
        The configuration object containing settings for the script.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.charmmff = kwargs.get('charmmff')
        assert self.charmmff is not None, 'charmmff must be provided to NAMDScripter'
        self.charmmff_config = kwargs.get('charmmff_config')
        self.charmrun = kwargs.get('charmrun')
        self.gpu_devices = kwargs.get('gpu_devices')
        self.local_ncpus = kwargs.get('local_ncpus')
        self.namd = kwargs.get('namd')
        self.namdgpu = kwargs.get('namd3gpu')
        self.namd_config = kwargs.get('namd_config')
        self.namd_type = kwargs.get('namd_type')
        self.namd_version = int(self.namd_config['namd-version'])
        logger.debug(f'Using NAMD v {self.namd_version}')
        self.ncpus = kwargs.get('ncpus')
        self.ngpus = kwargs.get('ngpus')
        self.slurmvars = kwargs.get('slurmvars')
        if self.namd_version == 2:
            self.namd_deprecates = {}
        else:
            self.namd_deprecates = kwargs.get('namd_deprecates')
        logger.debug(f'{self.ncpus} cpus are available for namd')
        self.default_ext = '.namd'

    def fetch_standard_charmm_parameters(self):
        """
        Fetch the standard CHARMM parameters from the configuration and copy them to the local directory.
        This method retrieves the parameters defined in the CHARMM force field configuration and copies them to
        the local directory for use in the NAMD script.
        
        Returns
        -------
        list
            A list of parameter files that have been copied to the local directory.
        """
        # logger.debug('Fetching standard CHARMMFF parameters')
        parameters_local = []
        for t in self.charmmff_config['standard']['prm'] + self.charmmff_config['standard']['str']:
            if not t in parameters_local:
                self.charmmff.copy_charmmfile_local(t)
                parameters_local.append(t)
        for t in self.charmmff_config['custom']['prm'] + self.charmmff_config['custom']['str']:
            if t not in parameters_local:
                self.charmmff.copy_charmmfile_local(t)
                parameters_local.append(t)
        # logger.debug(f'local parameters: {parameters_local}')
        return parameters_local

    def newscript(self, basename=None, addl_paramfiles=[]):
        """
        Initialize a new NAMD script with a specified basename and additional parameter files.
        If no basename is provided, a default script name is used.

        Parameters
        ----------
        basename : str, optional
            The base name for the script file. If not provided, a default name is used.
        addl_paramfiles : list, optional
            A list of additional parameter files to be included in the script. These files should not contain
            a path, but they should be recognizable .prm or .str files within the CHARMM force field or
            the custom CHARMM-format directories.
        """
        super().newscript(basename)
        self.scriptname = f'{basename}{self.default_ext}'
        self.banner('NAMD script')
        self.parameters = self.fetch_standard_charmm_parameters()
        for at in sorted(addl_paramfiles):
            if not at in self.parameters:
                self.charmmff.copy_charmmfile_local(at)
                self.parameters.append(at)
        for p in self.parameters:
            assert os.sep not in p
            self.addline(f'parameters {p}')

    def writescript(self, params, cpu_override=False):
        """
        Write the NAMD script based on the provided parameters.
        This method constructs the NAMD script by adding lines based on the parameters provided.
        If the NAMD type is 'gpu' and `cpu_override` is False, it will also add GPU-specific configurations.
        
        Parameters
        ----------
        params : dict
            A dictionary of parameters to include in the NAMD script.
        cpu_override : bool, optional
            If True, the script will be written with CPU-specific configurations, even if the NAMD type is 'gpu'.
        """
        # logger.debug(f'params: {params}')
        tailers = ['minimize', 'run', 'numsteps']
        for k, v in params.items():
            if k not in tailers:
                if type(v) == list:
                    for val in v:
                        if k == 'tcl':
                            self.addline(val)
                        else:
                            self.addline(f'{self.namd_deprecates.get(k, k)} {val}')
                else:
                    if k == 'tcl':
                        self.addline(v)
                    else:
                        self.addline(f'{self.namd_deprecates.get(k, k)} {v}')
        if self.namd_type == 'gpu' and not cpu_override:
            for k, v in self.namd_config['gpu-resident'].items():
                self.addline(f'{k} {v}')
        for t in tailers:
            if t in params:
                self.addline(f'{self.namd_deprecates.get(t, t)} {params[t]}')
        super().writescript()

    def runscript(self, **kwargs):
        """
        Run the NAMD script using the NAMD command line interface.
        This method constructs a command to execute NAMD with the specified script and options.

        Parameters
        ----------
        kwargs : dict
            A dictionary of options to be passed to the NAMD command. This can include:
            - ``local_execution_only``: If True, use the local CPU count instead of the configured CPU count.
            - ``single_gpu_only``: If True, use only one GPU device.
            - ``cpu_override``: If True, force the use of CPU settings even if the NAMD type is ``gpu``.
        """
        assert hasattr(self, 'scriptname'), f'No scriptname set.'
        if kwargs.get('local_execution_only', False):
            use_cpu_count = self.local_ncpus
        else:
            use_cpu_count = self.ncpus
        if kwargs.get('single_gpu_only', False):
            use_gpu_count = 1
            use_gpu_devices = '0'
        else:
            use_gpu_count = self.ngpus
            use_gpu_devices = self.gpu_devices
        if self.namd_type == 'cpu' or kwargs.get('cpu_override', False):
            c = Command(f'{self.charmrun} +p {use_cpu_count} {self.namd} {self.scriptname}')
        elif self.namd_type == 'gpu':
            if len(self.slurmvars) > 0:
                use_cpu_count = 8 if use_gpu_count == 1 else (use_gpu_count - 1) * 8 + 8 - (use_gpu_count - 1)
            else:
                use_cpu_count = self.local_ncpus
            c = Command(f'{self.namdgpu} +p{use_cpu_count} +setcpuaffinity +devices {use_gpu_devices} {self.scriptname}')
        self.logname = f'{self.basename}.log'
        self.logparser = NAMDLogParser(basename=self.basename)
        self.logparser.auxparser = NAMDxstParser(basename=f'{self.basename}')
        progress_struct = None
        if self.progress:
            progress_struct = NAMDProgress()
            self.logparser.enable_progress_bar(progress_struct)
        return c.run(logfile=self.logname, logparser=self.logparser)

