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
from .pipeline import PipelineContext
from .artifacts import TclScript, PDBFile, NAMDCoorFile, VMDLogFile, Artifact, ArtifactList

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
    index: int = 0
    taskname: str = 'BaseTask'
    pipeline: PipelineContext | None = None
    specs: dict = {}
    prior: BaseTask | None = None

    basename: str = ''
    subtaskcount: int = 0
    result: int = 0

    _extra_message: str = ''
    """
    An optional extra message that can be appended to the task's log messages.
    This can be used to provide additional context or information about the task's execution.
    """

    _init_msg_options: tuple[str] = ('INITIATED','STARTED','BEGUN','SET IN MOTION','KICKED OFF','LIT','SPANKED ON THE BOTTOM', 'LAUNCHED', 'KICKED OUT OF THE NEST', 'COMMENCED', 'ACTIVATED', 'UNLEASHED', 'ENGAGED', 'SPAWNED', 'BIRTHED', 'BORN', 'KICKED INTO GEAR', 'KICKED INTO ACTION', 'KICKED INTO HIGH GEAR', 'KICKED INTO OVERDRIVE', 'KICKED INTO HYPERDRIVE', 'KICKED INTO LIGHTSPEED')
    """
    A list of message options that are used to log the initiation of the task.
    """

    def __init__(self, index: int = 0, pipeline: PipelineContext = None, specs: dict = None, prior: BaseTask = None):
        """
        Constructor for the BaseTask class.
        """
        self.pipeline = pipeline if pipeline else PipelineContext()
        self.specs = specs if specs else {}
        self.prior = prior
        if prior:
            self.index = prior.index + 1
        else:
            self.index = index

        self.controller_index = self.pipeline.controller_index
        self.taskname = self.specs.get('taskname', f'{self.taskname}-{self.index:02d}')
        self.basename = ''
        self.subtaskcount = 0
        self.result = 0

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__} index={self.index} controller_index={self.controller_index} taskname={self.taskname}>"

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
        msg = 'initiated'
        if self._extra_message:
            msg += f' ({self._extra_message})'
        self.log_message(msg)
        self.result = self.do()
        if self.result == 0:
            msg = 'complete'
            if self._extra_message:
                msg += f' ({self._extra_message})'
            self.log_message(msg)
        else:
            logger.error(f'Task {self.taskname} failed with result {self.result}')

    def stash_current_artifact(self, key):
        """
        Stash the current artifact with the specified key into the history.
        This method moves the current artifact from the context to the history list.
        """
        self.pipeline.stash(key)

    def get_current_artifact_value(self, key):
        """ 
        Get the current artifact value with the specified key from the context.
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
        artifact=self.pipeline.get_artifact(key)
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
        artifact = self.pipeline.get_artifact(key)
        if artifact:
            return artifact.path if hasattr(artifact,'path') else artifact.value
        else:
            logger.debug(f'Artifact {key} not found')
        return default
    
    def get_current_artifact(self, key):
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
        return self.pipeline.get_artifact(key)

    def get_my_artifactfile_collection(self):
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
        all_my_artifacts, filelist_artifacts, file_artifacts = self.pipeline.get_artifact_collection_as_lists(produced_by=self)
        logger.debug(f'Found {len(all_my_artifacts)} artifacts, {len(filelist_artifacts)} file list artifacts, and {len(file_artifacts)} file artifacts for task {repr(self)}.')
        return file_artifacts

    def register_current_artifact(self, artifact: Artifact | ArtifactList, key: str = None):
        """
        Set the current artifact with the specified key in the context.

        Parameters
        ----------
        artifact : Artifact
            The artifact to register.

        """
        self.pipeline.register(artifact.stamp(self), key=key)

    def override_taskname(self, taskname):
        """
        Override the task name with a new name.  Certain situations may require a Controller to override the task name.
        This method is used to change the task name after the task has been created.

        Parameters
        ----------
        taskname : str
            The new task name to set.
        """
        logger.debug(f'Overriding task name {self.taskname} to {taskname}')
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
        extra = ''
        for k, v in kwargs.items():
            if v:
                extra += f' ({k}: {v})'
        mtoks = [x.strip() for x in [x.upper() for x in message.split()]]
        if not any([x in self._init_msg_options for x in mtoks]):
            extra += f' (result: {self.result})'
        logger.info(f'Controller {self.controller_index:02} Task {self.index:02} \'{self.taskname}\' {message} {extra}')

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

class VMDTask(BaseTask, ABC):
    """
    A base class for tasks that require VMD scripting.
    This class extends the BaseTask class and provides additional functionality for tasks that involve VMD scripting.
    It includes methods for converting coordinate files to PDB files, converting PDB files to coordinate files, and creating constraint PDB files.
    The VMDTask class is intended to be subclassed.
    """

    def coor_to_pdb(self):
        """
        Converts the coordinate file to a PDB file using VMD.
        This method creates a new VMD script to perform the conversion.
        It uses the ``namdbin2pdb`` Tcl proc to convert the coordinate file to a PDB file.
        The resulting PDB file is named based on the task's basename.
        """
        vm = self.pipeline.get_scripter('vmd')
        vm.newscript(f'{self.basename}-coor2pdb')
        psf = self.get_current_artifact_path('psf')
        if not psf:
            raise RuntimeError(f'No PSF file found for task {self.taskname}')
        coor = self.get_current_artifact_path('coor')
        if not coor:
            raise RuntimeError(f'No coordinate file found for task {self.taskname}')
        vm.addline(f'namdbin2pdb {psf.name} {coor.name} {self.basename}.pdb')
        vm.writescript()
        vm.runscript()
        self.register_current_artifact(TclScript(f'{self.basename}-coor2pdb'))
        self.register_current_artifact(VMDLogFile(f'{self.basename}-coor2pdb'))
        self.register_current_artifact(PDBFile(self.basename))

    def pdb_to_coor(self):
        """
        Converts the PDB file to a coordinate file using VMD.
        This method creates a new VMD script to perform the conversion.
        It uses the ``pdb2namdbin`` Tcl proc to convert the PDB file to a coordinate file.
        The resulting coordinate file is named based on the task's basename.
        """
        vm = self.pipeline.get_scripter('vmd')
        vm.newscript(f'{self.basename}-pdb2coor')
        pdb = self.get_current_artifact_path('pdb')
        if not pdb:
            raise RuntimeError(f'No PDB file found for task {self.taskname}')
        vm.addline(f'pdb2namdbin {pdb.name} {self.basename}.coor')
        vm.writescript()
        vm.runscript()
        self.register_current_artifact(TclScript(f'{self.basename}-pdb2coor'))
        self.register_current_artifact(VMDLogFile(f'{self.basename}-pdb2coor'))
        self.register_current_artifact(NAMDCoorFile(f'{self.basename}'))

    def make_constraint_pdb(self, specs, statekey='consref'):
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
        vm = self.pipeline.get_scripter('vmd')
        pdb = self.get_current_artifact_path('pdb')
        force_constant = specs.get('k', self.config['user']['namd']['harmonic']['spring_constant'])
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
        self.register_current_artifact(TclScript(f'{self.basename}-make-constraint-pdb'))
        vm.runscript()
        self.register_current_artifact(VMDLogFile(f'{self.basename}-make-constraint-pdb'))
        self.register_current_artifact(PDBFile(c_pdb), key=statekey)

