#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
A controller for the pestifer runtime.  Initialization of a Controller object generates
the configuration, from which the list of tasks is created.  The :meth:`do_tasks()` method executes
the tasks.
"""
import logging
import os

from .config import Config
from .pipeline import PipelineContext
from .basetask import BaseTask, TaskList

from ..scripters.filewriter import Filewriter
from ..scripters.namdscripter import  NAMDScripter
from ..scripters.packmolscripter import PackmolScripter
from ..scripters.psfgenscripter import PsfgenScripter
from ..scripters.tclscripters import VMDScripter

from ..tasks import task_classes, TerminateTask

logger = logging.getLogger(__name__)

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

    def __init__(self, config: Config, userspecs = {}, index = 0, terminate = True):
        self.index = index
        self.config = config
        self.config['user'].update(userspecs)

        logger.debug(f'task_classes {task_classes}')

        # set up the task list
        self.tasks = TaskList([])
        prior_task: BaseTask = None
        self.pipeline = PipelineContext(controller_index=self.index, 
                                        scripters={
                                            'psfgen': PsfgenScripter(self.config),
                                            'vmd': VMDScripter(self.config),
                                            'namd': NAMDScripter(self.config),
                                            'packmol': PackmolScripter(self.config),
                                            'data': Filewriter()
                                        }, 
                                        global_config=self.config['user'],
                                        resource_manager=self.config.RM)

        for idx, taskdict in enumerate(self.config['user'].get('tasks', [])):
            # Each task dictionary has a single keyword (the task name) and a value
            # that comprises the task specifications from the config
            assert len(taskdict) == 1, f"Task dictionary {taskdict} must have a single key-value pair"
            taskname = list(taskdict.keys())[0]
            task_specs = taskdict[taskname]
            task_specs['taskname'] = taskname
            logger.debug(f'{taskname}: {task_specs}')
            # specs={} if not config_specs else config_specs.copy()
            # Ensure the name of the task is among the implemented Tasks
            class_name = task_classes.get(taskname, None)
            # Create this Task instance:
            #   1. Get the class type
            #   2. Initialize the instance with 
            #      - specs: specifications from user config
            #      - taskname: the name of the task
            #      - self.config: the Controller configuration
            #      - self.scripters: the Controller's filewriters
            #      - prior_task: indentifier of prior task in task list
            Cls = task_classes[class_name]
            this_task = Cls(index = idx, pipeline = self.pipeline, specs = task_specs, prior = prior_task)
            # Append to the task list
            self.tasks.append(this_task)
            prior_task = this_task
        if len(self.tasks) > 0: logger.debug(f'last task currently is {self.tasks[-1].taskname}')
        # Add a "terminate" task by default if the user has not specified one
        if terminate and (len(self.tasks) == 0 or not isinstance(self.tasks[-1], TerminateTask)):
            # If the last task is not a TerminateTask, add one with default specs
            specs = self.config.make_default_specs('tasks','terminate')
            specs['taskname'] = 'terminate'
            logger.debug('Adding default terminate task')
            self.tasks.append(TerminateTask(index=len(self.tasks), pipeline=self.pipeline, specs=specs, prior=prior_task))
        ess='s' if len(self.tasks)>1 else ''
        logger.info(f'Running "{self.config["user"]["title"]}"')
        logger.info(f'Controller {self.index:02d} will execute {len(self.tasks)} task{ess}.')

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