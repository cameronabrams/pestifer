#Author: Cameron F. Abrams, <cfa22@drexel.edu>
"""A controller for the pestifer runtime
"""
import logging
logger=logging.getLogger(__name__)

from .config import Config
from .scriptwriters import Filewriter,Psfgen,VMD,NAMD2
from .tasks import *
from .util import *

class Controller:
    def __init__(self,userconfigfilename):
        # Read in the user configuration file and set up the overall Config
        self.config=Config(userfile=userconfigfilename)

        # Set up the file writers
        self.writers={
            'psfgen': Psfgen(self.config),
            'vmd':    VMD(self.config),
            'namd2':  NAMD2(self.config),
            'data':   Filewriter()
        }

        # set up the task list
        task_classes=inspect_classes('pestifer.tasks')
        self.tasks=[]
        prior_task=None
        for taskdict in self.config['user']['tasks']:
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
        logger.debug(f'last task currently is {self.tasks[-1].taskname}')
        # Add a "terminate" task by default if the user has not specified one
        if len(self.tasks)==0 or not self.tasks[-1].taskname=='terminate':
            specs=self.config.make_default_specs('tasks','terminate')
            logger.debug('Adding default terminate task')
            self.tasks.append(TerminateTask(specs,'terminate',self.config['base'],self.writers,prior_task))
        logger.info(f'Controller will execute {len(self.tasks)} task(s).')

    def do_tasks(self):
        # Execute each task in series
        for task in self.tasks:
            task.do()

    def write_complete_config(self,filename='00-complete-user.yaml'):
        self.config.dump_user(filename=filename)