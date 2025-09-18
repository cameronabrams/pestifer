#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
A controller for the pestifer runtime.  Initialization of a Controller object generates
the configuration, from which the list of tasks is created.  The :meth:`do_tasks()` method executes
the tasks.
"""
import logging

from dataclasses import dataclass

from pestifer.tasks.validate import ValidateTask

from .config import Config
from .pipeline import PipelineContext
from ..tasks.taskcollections import TaskList
from ..tasks import TerminateTask
from ..util.stringthings import plu, my_logger
from ..util.util import running_under_pytest
from ..util._goldenmode import report_example_id

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

        # if report_example_id() and running_under_pytest() and terminate:
        #     # truncate task list after latest ValidateTask, if one exists
        #     # remove the TerminateTask and put on a custom one so that 
        #     # desired outputs are retained
        #     from ..tasks.pytest_buildterminate import Pytest_BuildTerminateTask
        #     validate_task_idx = None
        #     for i, task in enumerate(self.tasks):
        #         if isinstance(task, ValidateTask):
        #             validate_task_idx = i
        #     if validate_task_idx is not None:
        #         self.tasks = self.tasks[:validate_task_idx+1]
        #     else:
        #         self.tasks = self.tasks[:-1] # remove the original terminate task
        #     # add the special pytest terminate task
        #     self.tasks.append(Pytest_BuildTerminateTask(specs={'prior': self.tasks[-1]}, index=len(self.tasks)))

        self.pipeline = PipelineContext(controller_index=self.index)
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
        logger.debug(f'Controller {self.index:02d} provisioning {len(self.tasks)} task{plu(len(self.tasks))} with packet:')
        my_logger(packet, logger.debug)
        for task in self.tasks:
            logger.debug('*'*70)
            logger.debug(f'Task {self.index:02d}:{task.index:02d} ** {task.taskname:>20s} **')
            task.provision(packet)
            if task.specs.get('requires_subcontroller', False):
                logger.debug('*'*70)
                logger.debug(f'Controller {self.index:02d} spawning subcontroller for Task {task.index:02d} ** {task.taskname:>20s} **')
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
        task_durations = 0
        for task in self.tasks:
            returned_result = task.execute()
            task_report[task.index] = dict(taskname=task.taskname, taskindex=task.index, result=returned_result)
            task_durations += task.duration
            if task.result != 0:
                logger.warning(f'Task {task.taskname} failed; task.result {task.result} returned result {returned_result} controller is aborted.')
                break
        for task in self.tasks:
            task_report[task.index]['duration'] = task.duration
            task_report[task.index]['duration_frac'] = task.duration/task_durations if task_durations > 0 else 0
        return task_report

    def write_complete_config(self, filename='complete-user.yaml'):
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