# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
This module defines the :class:`BaseTask` class, which serves as a base class for all task types in the pestifer framework.
The :class:`BaseTask` class provides a common interface and some basic functionality for tasks, including task initialization,
logging, state management, and file handling. It is intended to be subclassed by specific task types that implement their own
specific behavior and functionality.
The :class:`BaseTask` class includes methods for task execution, state variable management, file collection, and task name management.
It also provides a mechanism for inheriting state from prior tasks, saving the task state to a YAML file, and copying the state to files based on specified extensions.
The class also includes methods for converting coordinate files to PDB files and vice versa, as well as creating constraint PDB files based on task specifications.
The :class:`BaseTask` class is designed to be flexible and extensible, allowing for the creation of various task types that can perform different operations on molecular structures.

Available tasks that inherit from :class:`BaseTask` include all those in the :mod:`pestifer.tasks` subpackage.

"""
from __future__ import annotations

import logging

from abc import ABC, abstractmethod
from time import perf_counter
from typing import TYPE_CHECKING

from ..core.artifacts import *
from ..core.command import Command
from ..core.resourcemanager import ResourceManager
from ..core.pipeline import PipelineContext

from ..scripters import GenericScripter, VMDScripter
from ..util.util import hmsf

if TYPE_CHECKING:
    from ..core.controller import Controller

logger = logging.getLogger(__name__)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

class BaseTask(ABC):
    """ 
    A base class for Tasks. This class is intended to be subclassed by specific task types.
    It provides a common interface and some basic functionality for tasks.  All tasks are
    assigned to a PipelineContext owned by a Controller.

    Parameters
    ----------
    config_specs : dict, optional
        A dictionary of configuration specifications. This can include task-specific settings.
    controller_specs : dict, optional
        A dictionary of specifications provided by the controller that created this task.
        It can include attributes like `prior`, `index`, `scripters`, `taskname`, and `controller_index`.

    """

    _yaml_header = 'task'  # used to identify the task in YAML files

    _init_msg_options: tuple[str] = ('INITIATED', 'STARTED', 'BEGUN', 'SET IN MOTION', 'KICKED OFF', 'LIT', \
                                     'SPANKED ON THE BOTTOM', 'LAUNCHED', 'KICKED OUT OF THE NEST', 'COMMENCED', \
                                     'ACTIVATED', 'UNLEASHED', 'ENGAGED', 'SPAWNED', 'BIRTHED', 'BORN', \
                                     'KICKED INTO GEAR', 'KICKED INTO ACTION', 'KICKED INTO HIGH GEAR', \
                                     'KICKED INTO OVERDRIVE', 'KICKED INTO HYPERDRIVE', 'KICKED INTO LIGHTSPEED')
    """
    A list of message options that are used to log the initiation of the task.
    """

    def __init__(self, specs: dict = None, **kwargs):
        """
        Constructor for the BaseTask class.
        """
        self.specs: dict = specs if specs else {}
        self.taskname: str = kwargs.get('taskname', f'{self._yaml_header}')
        self.index: int = kwargs.get('index', None)
        self.provisions: dict = None
        self.subcontroller: 'Controller' = None
        self.prior: BaseTask = None

        self.basename: str = ''
        self.subtaskcount: int = 0
        self.result: int = 0
        self.extra_message: str = ''

        self.pytest_skip_after_to_terminate = self.specs.get('pytest_skip_after_to_terminate', False)

    def provision(self, packet: dict = None):
        """
        Provision the task with necessary resources.
        This method is called to set up the task with the required resources and context.
        It can be overridden by subclasses to provide additional provisioning logic.

        Parameters
        ----------
        packet : dict, optional
            A dictionary of provisions to be used for the task. If not provided, the existing provisions will be used.
        """
        if packet:
            if self.provisions is None:
                self.provisions = {}
            self.provisions.update(packet)
        
        # set some shortcuts for commonly used provisions
        self.resource_manager: ResourceManager = self.provisions.get('resource_manager', None)
        self.scripters: dict[GenericScripter] = self.provisions.get('scripters', {})
        self.controller_index: int = self.provisions.get('controller_index', 0)
        self.pipeline: PipelineContext = self.provisions.get('pipeline', None)

    @property
    def is_provisioned(self) -> bool:
        """
        Check if the task has been provisioned with necessary resources.
        
        Returns
        -------
        bool
            True if the task has been provisioned, False otherwise.
        """
        return bool(self.provisions)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} index={self.index} controller_index={self.controller_index} taskname={self.taskname}>"

    def get_scripter(self, name: str):
        """
        Get a scripter by name from the pipeline's scripters.
        
        Parameters
        ----------
        name : str
            The name of the scripter to retrieve.

        Returns
        -------
        Scripter
            The scripter instance associated with the given name.
        """
        query_result = self.scripters.get(name, None)
        if query_result is None:
            logger.debug(f'Scripter {name} not found')
            if not self.is_provisioned:
                logger.warning(f'This task has not been provisioned.')
        return query_result

    @abstractmethod
    def do(self) -> int:
        """
        This a stub method that should be overridden by subclasses.
        It is intended to be the main method that performs the task's operations.
        Subclasses should implement this method to define the specific behavior of the task.
        This method is called when the task is executed.
        """
        return -1

    def execute(self) -> int:
        """
        Execute the task.
        This method calls the `do` method, which should be implemented by subclasses to perform the task's operations.
        It also logs the initiation and completion of the task.
        """
        if not self.is_provisioned:
            logger.warning(f'Task {self.taskname} is not provisioned.')
            return -1
        msg = 'initiated'
        if self.extra_message:
            msg += f' ({self.extra_message})'
        self.log_message(msg)
        t1 = perf_counter()
        self.result = self.do()
        t2 = perf_counter()
        if self.result == 0:
            msg = 'completed'
        else:
            msg = 'failed'
        if self.extra_message:
            msg += f' ({self.extra_message})'
        self.duration = t2 - t1
        msg += f' in {hmsf(self.duration)}'
        self.log_message(msg)
        return self.result
        
    def stash_current_artifact(self, key):
        """
        Stash the current artifact with the specified key into the history.
        This method moves the current artifact from the context to the history list.
        """
        self.pipeline.stash(key)

    def get_current_artifact_data(self, key):
        """ 
        Get the current artifact data with the specified key from the context.
        This method retrieves the artifact from the context that matches the specified key.
        
        Parameters
        ----------
        key : str
            The key of the artifact to retrieve from the context.

        Returns
        -------
        Artifact
            The artifact associated with the specified key, or None if not found.
        """
        artifact: Artifact | FileArtifact = self.pipeline.get_current_artifact(key)
        if artifact:
            return artifact.data
        else:
            logger.debug(f'Artifact {key} not found')
        return None

    def get_current_artifact_path(self, key, default=None):
        """
        Get the current artifact path with the specified key from the context.
        This method retrieves the artifact from the context that matches the specified key.
        
        Parameters
        ----------
        key : str
            The key of the artifact to retrieve from the context.

        Returns
        -------
        Artifact
            The artifact associated with the specified key, or None if not found.
        """
        artifact: Artifact | FileArtifact = self.pipeline.get_current_artifact(key)
        if artifact:
            return artifact.path if hasattr(artifact,'path') else artifact.data
        else:
            logger.debug(f'Artifact {key} not found')
        return default
    
    def get_current_artifact(self, key, **kwargs):
        """
        Get the current artifact with the specified key from the context.
        This method retrieves the artifact from the context that matches the specified key.

        Parameters
        ----------
        key : str
            The key of the artifact to retrieve from the context.

        Returns
        -------
        Artifact
            The artifact associated with the specified key, or None if not found.
        """
        return self.pipeline.get_current_artifact(key, **kwargs)

    def get_my_artifactfile_collection(self) -> FileArtifactList:
        """
        Get a collection of artifact files produced by the current task.
        This method retrieves all artifact files that have been registered in the context
        and are associated with the current task.

        Returns
        -------
        list
            A list of artifact files associated with the current task.
        """
        logger.debug(f'Getting artifact file collection for task {repr(self)}')
        all_my_artifact_files = FileArtifactList(self.pipeline.get_all_file_artifacts(produced_by=self))
        logger.debug(f'Found {len(list(all_my_artifact_files))} artifacts for task {repr(self)}.')
        return all_my_artifact_files

    def register(self, data: object, key: str, artifact_type = DataArtifact, **kwargs):
        """
        Register data as an artifact in the pipeline context.

        Parameters
        ----------
        data : object
            The data to be registered as an artifact.
        key : str
            The key under which the artifact will be registered.
        **kwargs : dict
            Additional keyword arguments that may include metadata or other information
        """
        return self.pipeline.register(data, key, self, artifact_type=artifact_type, **kwargs)

    def register_if_exists(self, data: object, requestor: object, artifact_type: type = FileArtifact, **kwargs) -> FileArtifact | None:
        return self.pipeline.register_if_exists(data, requestor, artifact_type=artifact_type, **kwargs)
    
    def override_taskname(self, taskname):
        """
        Override the task name with a new name.  Certain situations may require a Controller to override the task name.
        This method is used to change the task name after the task has been created.

        Parameters
        ----------
        taskname : str
            The new task name to set.
        """
        # logger.debug(f'Overriding task name {self.taskname} to {taskname}')
        self.taskname = taskname

    def __str__(self):
        return repr(self)

    def log_message(self, message, **kwargs):
        """ 
        Logs a message with the task's index and name.

        Parameters
        ----------
        message : str
            The message to log.
        **kwargs : dict
            Additional keyword arguments that can be used to add extra information to the log message.
            These are typically used to provide additional context or details about the task's execution.
        """
        if self.index is None:
            logger.info(f'Task is unindexed: {message}')
            return
        extra = ''
        for k, v in kwargs.items():
            if v:
                extra += f' ({k}: {v})'
        mtoks = [x.strip() for x in [x.upper() for x in message.split()]]
        if not any([x in self._init_msg_options for x in mtoks]):
            extra += f' (result: {self.result})'
        logger.info(f'Controller {self.controller_index:02d} Task {self.index:02d} \'{self.taskname}\' {message} {extra}')

    def next_basename(self, extra_label: str = ''):
        """
        Generates a new basename for the task based on the controller index, task index, subtask count, and task name.
        If the ``specs`` dictionary contains a 'basename' key, it overrides the default basename.
        The basename is constructed in the format: ``(controller_index)-(task_index)-(subtask_count)_(taskname)[-(label)]``.
        If ``extra_label`` is provided, it is used as the label.  For example, next_basename('mylabel') would result in a basename like ``01-02-03_taskname-mylabel``.
        """

        label = ''
        if extra_label:
            label = f'-{extra_label}'
        default_basename = f'{self.controller_index:02d}-{self.index:02d}-{self.subtaskcount:02d}_{self.taskname}{label}'
        overwrite_basename = self.specs.get('basename', None)  # user put a basename in the task specs in the YAML input
        basename = overwrite_basename if overwrite_basename else default_basename
        self.basename = basename
        self.subtaskcount += 1

    def import_artifacts(self, pipeline: PipelineContext):
        """
        Imports artifacts from the pipeline context into the task.
        """
        self.pipeline.import_artifacts(pipeline)

class VMDTask(BaseTask, ABC):
    """
    A base class for tasks that require VMD scripting.
    This class extends the BaseTask class and provides additional functionality for tasks that involve VMD scripting.
    It includes methods for converting coordinate files to PDB files, converting PDB files to coordinate files, and creating constraint PDB files. The VMDTask class is intended to be subclassed.
    """

    def coor_to_pdb(self, coorfilename: str, psffilename: str) -> str:
        """
        Converts a namdbin coordinate file to a PDB file using catdcd.
        """
        if not Path(coorfilename).exists():
            raise FileNotFoundError(f'Coordinate file {coorfilename} not found')
        c = Command(f'catdcd -o {self.basename}.pdb -otype pdb -s {psffilename} -stype psf -namdbin {coorfilename}')
        c.run()
        with open(f'{self.basename}.pdb', 'r') as f:
            pdb_lines = f.readlines()
        with open(f'{self.basename}.pdb', 'w') as f:
            for line in pdb_lines:
                if not line.startswith('REMARK') and not line.startswith('CRYST'):
                    f.write(line)
        return f'{self.basename}.pdb'

    def pdb_to_coor(self, pdbfilename: str) -> str:
        """
        Converts a PDB file to a namdbin coordinate file using catdcd.
        """
        c = Command(f'catdcd -o {self.basename}.coor -otype namdbin -stype psf -pdb {pdbfilename}')
        c.run()
        return f'{self.basename}.coor'

    def make_constraint_pdb(self, specs: dict, statekey: str = 'consref'):
        """
        Creates a PDB file with constraints based on the specifications provided.
        This method generates a VMD script that selects atoms based on the specifications and writes a PDB file with the specified constraints.
        The resulting PDB file is named based on the task's basename and the state key.
        This allows the task to keep track of the constraints PDB file in its state.
        The method uses the VMD script writer to create and run the script that generates the constraints PDB file.
        The script selects all atoms, sets their occupancy to 0.0, selects the constrained atoms based on the specifications, and sets their occupancy to the specified force constant.
        Finally, it writes the PDB file with the specified name.
        If the 'consref' key is not provided in the specifications, it defaults to ``f'{self.basename}-constraints.pdb'``.
        The resulting PDB file is then stored in the task's state variables under the specified statekey.
        If the specified atoms are not found in the PDB, an error will be raised.

        Parameters
        ----------
        specs : dict
            A dictionary containing the specifications for the constraints. It should include:

                - ``k``: The force constant for the constraints (default is taken from the config).
                - ``atoms``: A string defining the atoms to be constrained (default is ``all``).
                - ``consref``: The name of the output PDB file for constraints (default is ``f'{self.basename}-constraints.pdb'``).

        statekey : str, optional
            The key under which the resulting PDB file will be stored in the state variables. Default is ``consref``.
        """
        vm: VMDScripter = self.get_scripter('vmd')
        state: StateArtifacts = self.get_current_artifact('state')
        pdb: Path = state.pdb
        # 'namd_global_config': self.config['user']['namd'],
        force_constant = specs.get('k', self.namd_global_config['harmonic']['spring_constant'])
        constrained_atoms_def = specs.get('atoms', 'all')
        logger.debug(f'constraint spec: {specs["atoms"]}')
        c_pdb = specs.get('consref', '')
        if not c_pdb:
            c_pdb = f'{self.basename}-constraints.pdb'
        vm.newscript(f'{self.basename}-make-constraint-pdb')
        vm.addline(f'mol new {pdb.name}')
        vm.addline(f'set a [atomselect top all]')
        vm.addline(f'$a set occupancy 0.0')
        vm.addline(f'set c [atomselect top "{constrained_atoms_def}"]')
        vm.addline(f'$c set occupancy {force_constant}')
        vm.addline(f'$a writepdb {c_pdb}')
        vm.writescript()
        self.register(f'{self.basename}-make-constraint-pdb', key='tcl', artifact_type=TclScriptArtifact)
        vm.runscript()
        self.register(f'{self.basename}-make-constraint-pdb', key='log', artifact_type=VMDLogFileArtifact)
        self.register(c_pdb, key=statekey, artifact_type=PDBFileArtifact)

