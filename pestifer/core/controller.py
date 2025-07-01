#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" 
A controller for the pestifer runtime.  Initialization of a Controller object generates
the configuration, from which the list of tasks is created.  The :meth:`do_tasks()` method executes
the tasks.
"""
import logging
import os

from .config import Config
from .scripters import Filewriter, PsfgenScripter, VMDScripter, NAMDScripter
from ..tasks.terminate import TerminateTask
from ..util.util import inspect_package_dir

from .. import tasks

logger=logging.getLogger(__name__)

class Controller:
    """ 
    Controller class for managing the execution of tasks in the Pestifer runtime.
    This class initializes with a configuration and a list of user-specified tasks.
    It sets up the necessary file writers and creates a list of tasks to be executed.
    The tasks are executed in the order they are defined, and the results of each task
    are collected in a report. If the last task is not a :class:`pestifer.tasks.terminate.TerminateTask` task,
    a default :class:`pestifer.tasks.terminate.TerminateTask` is added to ensure proper cleanup.

    Parameters
    ----------
    config : Config
        The configuration object containing user specifications and settings.
    userspecs : dict, optional
        A dictionary of user-specific configurations that will update the base configuration.
    index : int, optional
        An index for the controller instance, useful for identifying multiple controllers in a run.
    """
    def __init__(self,config:Config,userspecs={},index=0):
        self.index=index
        self.config=config
        if userspecs:
            self.config['user'].update(userspecs)

        # Set up the scripters
        self.scripters={
            'psfgen': PsfgenScripter(self.config),
            'vmd':    VMDScripter(self.config),
            'namd':   NAMDScripter(self.config),
            'data':   Filewriter()
        }

        task_classes,dum=inspect_package_dir(os.path.dirname(tasks.__file__))
        logger.debug(f'task_classes {task_classes}')

        # set up the task list
        self.tasks=[]
        prior_task=None
        for taskdict in self.config['user'].get('tasks',[]):
            # Each task dictionary has a single keyword (the task name) and a value
            # that comprises the task specifications from the config
            assert len(taskdict)==1, f"Task dictionary {taskdict} must have a single key-value pair"
            taskname=list(taskdict.keys())[0]
            config_specs=taskdict[taskname]
            logger.debug(f'{taskname}: {config_specs}')
            # specs={} if not config_specs else config_specs.copy()
            # Ensure the name of the task is among the implemented Tasks
            class_name=[name for name,cls in task_classes.items() if cls.yaml_header==taskname][0]
            # Create this Task instance:
            #   1. Get the class type
            #   2. Initialize the instance with 
            #      - specs: specifications from user config
            #      - taskname: the name of the task
            #      - self.config: the Controller configuration
            #      - self.scripters: the Controller's filewriters
            #      - prior_task: indentifier of prior task in task list
            Cls=task_classes[class_name]
            this_task=Cls(config_specs,dict(controller_index=self.index,taskname=taskname,config=self.config,scripters=self.scripters,prior=prior_task))
            # Append to the task list
            self.tasks.append(this_task)
            prior_task=this_task
        if len(self.tasks)>0: logger.debug(f'last task currently is {self.tasks[-1].taskname}')
        # Add a "terminate" task by default if the user has not specified one
        if len(self.tasks)==0 or not self.tasks[-1].taskname=='terminate':
            specs=self.config.make_default_specs('tasks','terminate')
            logger.debug('Adding default terminate task')
            self.tasks.append(TerminateTask(specs,dict(
                controller_index=self.index,
                taskname='terminate',
                config=self.config['base'],
                scripters=self.scripters,
                prior=prior_task)))
        ess='s' if len(self.tasks)>1 else ''
        logger.info(f'Run title: "{self.config["user"]["title"]}"')
        logger.info(f'Controller {self.index:02d} will execute {len(self.tasks)} task{ess}.')

    def do_tasks(self):
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
        task_report={}
        for task in self.tasks:
            returned_result=task.do()
            task_report[task.index]=dict(taskname=task.taskname,taskindex=task.index,result=returned_result)
            if task.result!=0:
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