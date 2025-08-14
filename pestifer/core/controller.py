#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
A controller for the pestifer runtime.  Initialization of a Controller object generates
the configuration, from which the list of tasks is created.  The :meth:`do_tasks()` method executes
the tasks.
"""
import logging

from dataclasses import dataclass

from .config import Config
from .pipeline import PipelineContext
from ..tasks.taskcollections import TaskList
from ..tasks import TerminateTask
from ..util.stringthings import plu

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

    def configure(self, config: Config, userspecs: dict = {}, index: int = 0, terminate: bool = True):
        self.index = index
        self.config = config
        self.config['user'].update(userspecs)

        # set up the task list
        self.tasks = TaskList.from_yaml(self.config['user'].get('tasks', []))
        if terminate and (len(self.tasks) == 0 or not isinstance(self.tasks[-1], TerminateTask)):
            # If the last task is not a TerminateTask, add one with default specs
            specs = self.config.make_default_specs('tasks','terminate')
            specs['taskname'] = 'terminate'
            specs['index'] = len(self.tasks)
            specs['cleanup'] = False
            logger.debug('Adding default terminate task')
            self.tasks.append(TerminateTask(specs=specs))

        self.pipeline = PipelineContext(controller_index = self.index)
        self.provision_tasks()

        logger.info(f'Running "{self.config["user"]["title"]}"')
        logger.info(f'Controller {self.index:02d} will execute {len(self.tasks)} task{plu(len(self.tasks))}.')
        return self

    @classmethod
    def spawn_subcontroller(cls, progenitor: 'Controller') -> 'Controller':
        subcontroller = cls()
        subcontroller.config = progenitor.config
        subcontroller.index = progenitor.index + 1
        return subcontroller
    
    def reconfigure_tasks(self, tasks: list[dict]):
        self.tasks = TaskList.from_yaml(tasks)
        self.provision_tasks()
        logger.info(f'Reconfigured tasks in "{self.config["user"]["title"]}"')
        logger.info(f'Controller {self.index:02d} will execute {len(self.tasks)} task{plu(len(self.tasks))}.')

    def provision_tasks(self):
        packet = {
            'controller_index': self.index,
            'subcontroller': Controller.spawn_subcontroller(self),
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
            task.provision(packet)

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