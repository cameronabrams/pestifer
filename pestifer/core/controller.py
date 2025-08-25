#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
A controller for the pestifer runtime.  Initialization of a Controller object generates
the configuration, from which the list of tasks is created.  The :meth:`do_tasks()` method executes
the tasks.
"""
import logging

from dataclasses import dataclass

from pestifer.tasks.make_membrane_system import MakeMembraneSystemTask
from pestifer.tasks.mdtask import MDTask
from pestifer.tasks.psfgen import PsfgenTask
from pestifer.tasks.validate import ValidateTask

from .config import Config
from .pipeline import PipelineContext
from ..tasks.taskcollections import TaskList
from ..tasks import TerminateTask
from ..util.stringthings import plu
from ..util.util import running_under_pytest

logger = logging.getLogger(__name__)

@dataclass
class Controller:
    """ 
    Controller class for managing the execution of tasks in the Pestifer runtime.
    This class initializes with a configuration and a list of user-specified tasks.
    It sets up the necessary file writers and creates a list of tasks to be executed.
    The tasks are executed in the order they are defined, and the results of each task
    are collected in a report. If the last task is not a :class:`pestifer.tasks.terminate.TerminateTask` task,
    a default :class:`pestifer.tasks.terminate.TerminateTask` is added to ensure proper termination.

    Parameters
    ----------
    config : Config
        The configuration object containing user specifications and settings.
    userspecs : dict, optional
        A dictionary of user-specific configurations provided at run-time that will update the base configuration.
    index : int, optional
        An index for the controller instance, useful for identifying multiple controllers in a run.
    terminate : bool, optional
        If True, a default terminate task will be added if the last task is not already a TerminateTask.
        This ensures that the runtime will always have a way to clean up and terminate properly.
        Default is True.
    """

    def __init__(self):
        self.index = 0
        self.config = None
        self.tasks = None
        self.pipeline = None
        self.parent = None

    def configure(self, config: Config, userspecs: dict = {}, index: int = 0, terminate: bool = True):
        self.index = index
        self.config = config
        self.config.update_user(userspecs)

        # set up the task list
        self.tasks = TaskList.from_yaml(self.config['user'].get('tasks', []))
        if terminate and (len(self.tasks) == 0 or not isinstance(self.tasks[-1], TerminateTask)):
            # If the last task is not a TerminateTask, add one with default specs
            specs = self.config.make_default_specs('tasks','terminate')
            logger.debug('Adding default terminate task')
            self.tasks.append(TerminateTask(specs=specs, index=len(self.tasks)))

        # if this is being run by pytest, find the first task in the list for which
        # the pytest_skip_after_to_terminate is true, and delete all tasks between it and terminate.
        if running_under_pytest():
            psfgen_idx = None
            make_membrane_idx = None
            for i, task in enumerate(self.tasks):
                if isinstance(task, PsfgenTask):
                    psfgen_idx = i
                elif isinstance(task, MakeMembraneSystemTask):
                    make_membrane_idx = i
            if make_membrane_idx is not None:
                task_idx = make_membrane_idx + 1
                if isinstance(self.tasks[task_idx], ValidateTask):
                    task_idx += 1
                if isinstance(self.tasks[task_idx], MDTask):
                    task_idx += 1
                self.tasks = self.tasks[:task_idx] + [self.tasks[-1]]
                self.tasks[-1].prior = self.tasks[task_idx-1]
            elif psfgen_idx is not None:
                task_idx = psfgen_idx + 1
                if isinstance(self.tasks[task_idx], ValidateTask):
                    task_idx += 1
                if isinstance(self.tasks[task_idx], MDTask):
                    task_idx += 1
                self.tasks = self.tasks[:task_idx] + [self.tasks[-1]]
                self.tasks[-1].prior = self.tasks[task_idx-1]

        self.pipeline = PipelineContext(controller_index = self.index)
        self.provision_tasks()

        logger.info(f'Running "{self.config["user"]["title"]}"')
        logger.info(f'Controller {self.index:02d} will execute {len(self.tasks)} task{plu(len(self.tasks))}.')
        return self

    @classmethod
    def spawn_subcontroller(cls, progenitor: 'Controller') -> 'Controller':
        logger.debug(f'Spawning subcontroller from {progenitor.index:02d}')
        subcontroller = cls()
        subcontroller.configure(progenitor.config.taskless_subconfig(), index=progenitor.index+1, terminate=False)
        subcontroller.parent = progenitor
        return subcontroller
    
    def reconfigure_tasks(self, tasks: list[dict]):
        new_data = dict(tasks=tasks)
        self.config.update_user(new_data)
        self.tasks = TaskList.from_yaml(self.config['user'].get('tasks', []))
        self.provision_tasks()
        logger.info(f'Reconfigured tasks in "{self.config["user"]["title"]}"')
        logger.info(f'Controller {self.index:02d} will execute {len(self.tasks)} task{plu(len(self.tasks))}.')

    def provision_tasks(self):
        packet = {
            'controller_index': self.index,
            'pipeline': self.pipeline,
            'resource_manager': self.config.RM,
            'scripters': self.config.scripters,
            'processor-type': self.config['user']['namd']['processor-type'],
            'shell-commands': self.config.shell_commands,
            'progress-flag': self.config.use_terminal_progress,
            'namd_global_config': self.config['user']['namd'],
        }
        logger.debug(f'Provisioning tasks with packet: {packet}')
        for task in self.tasks:
            logger.debug(f'  -> task {task.index} {task.taskname}')
            task.provision(packet)
            if task.specs.get('requires_subcontroller', False):
                logger.debug(f'    -> spawning subcontroller for task {task.index} {task.taskname}')
                task.subcontroller = Controller.spawn_subcontroller(self)

    def do_tasks(self) -> dict:
        """
        Execute the tasks in the order they were defined.

        This method iterates through the list of tasks, executing each one in turn.
        It collects the results of each task in a report dictionary, which maps task indices to their
        names, indices, and results. If any task fails (i.e., returns a non-zero result),
        a warning is logged, and the execution of subsequent tasks is aborted.

        Returns
        -------
        dict
            A dictionary containing the results of each task, with task indices as keys.
            Each value is itself a dictionary with the following keys:

            - ``taskname``: The name of the task
            - ``taskindex``: The index of the task
            - ``result``: The task's return code
        """
        task_report = {}
        for task in self.tasks:
            returned_result = task.execute()
            task_report[task.index] = dict(taskname=task.taskname, taskindex=task.index, result=returned_result)
            if task.result != 0:
                logger.warning(f'Task {task.taskname} failed; task.result {task.result} returned result {returned_result} controller is aborted.')
                break
        return task_report

    def write_complete_config(self,filename='complete-user.yaml'):
        """ 
        Write the complete user configuration to a YAML file.
        This method dumps the user configuration, including all modifications made during the execution of tasks,
        to a specified YAML file. The default filename is 'complete-user.yaml'.
        
        Parameters
        ----------
        filename : str, optional
            The name of the file to which the user configuration will be written. Default is 'complete-user.yaml'.
        """
        self.config.dump_user(filename=filename)