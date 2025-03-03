#Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" A controller for the pestifer runtime.  Initialization of a Controller object generates
the configuration, from which the list of tasks is created.  The do_tasks() method executes
the tasks.
"""
import logging
import os

from .basetask import BaseTask
from .config import Config
from .scriptwriters import Filewriter,Psfgen,VMD,NAMD
from .tasks.terminate import TerminateTask
from .util.util import inspect_package_dir

from . import tasks

logger=logging.getLogger(__name__)

class Controller:
    """ A class for controlling the execution of `Tasks`.
    
    """
    def __init__(self,userconfigfilename='',userspecs={},**kwargs):
        # Read in the user configuration file and set up the overall Config
        self.config=Config(userfile=userconfigfilename,**kwargs)
        if userspecs:
            self.config['user'].update(userspecs)

        # Set up the file writers
        self.writers={
            'psfgen': Psfgen(self.config),
            'vmd':    VMD(self.config),
            'namd':   NAMD(self.config),
            'data':   Filewriter()
        }

        task_classes,dum=inspect_package_dir(os.path.dirname(tasks.__file__))
        logger.debug(f'task_classes {task_classes}')

        # set up the task list
        self.tasks=[]
        prior_task=None
        BaseTask._taskcount=0
        for taskdict in self.config['user'].get('tasks',[]):
            # Each task dictionary has a single keyword (the task name) and a value
            # that comprises the task specifications
            assert len(taskdict)==1
            taskname=list(taskdict.keys())[0]
            specval=taskdict[taskname]
            logger.debug(f'{taskname}: {specval}')
            specs={} if not specval else specval.copy()
            # Ensure the name of the task is among the implemented Tasks
            class_name=[name for name,cls in task_classes.items() if cls.yaml_header==taskname][0]
            # Create this Task instance:
            #   1. Get the class type
            #   2. Initialize the instance with 
            #      - specs: specifications from user config
            #      - taskname: the name of the task
            #      - self.config: the Controller configuration
            #      - self.writers: the Controller's filewriters
            #      - prior_task: indentifier of prior task in task list
            Cls=task_classes[class_name]
            this_task=Cls(specs,taskname,self.config,self.writers,prior_task)
            # Append to the task list
            self.tasks.append(this_task)
            prior_task=this_task
        if len(self.tasks)>0: logger.debug(f'last task currently is {self.tasks[-1].taskname}')
        # Add a "terminate" task by default if the user has not specified one
        if len(self.tasks)==0 or not self.tasks[-1].taskname=='terminate':
            specs=self.config.make_default_specs('tasks','terminate')
            logger.debug('Adding default terminate task')
            self.tasks.append(TerminateTask(specs,'terminate',self.config['base'],self.writers,prior_task))
        ess='s' if len(self.tasks)>1 else ''
        logger.info(f'Run title: "{self.config["user"]["title"]}"')
        logger.info(f'Controller will execute {len(self.tasks)} task{ess}.')

    def do_tasks(self):
        # Execute each task in series
        task_report={}
        for task in self.tasks:
            returned_result=task.do()
            task_report[task.index]=dict(taskname=task.taskname,taskindex=task.index,result=returned_result)
            if task.result!=0:
                logger.warning(f'Task {task.taskname} failed; task.result {task.result} returned result {returned_result} controller is aborted.')
                break
        return task_report

    def write_complete_config(self,filename='00-complete-user.yaml'):
        self.config.dump_user(filename=filename)