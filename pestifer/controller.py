"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

from .config import Config
from .scriptwriters import Filewriter,Psfgen,VMD,NAMD2
from .tasks import *
from .util import *

class Controller:
    def __init__(self,userconfigfilename):
        self.config=Config(userconfigfilename)
        logger.debug(f'New controller: rcsb format {self.config["rcsb_file_format"]}')
        self.writers={
            'psfgen': Psfgen(self.config),
            'vmd':    VMD(self.config),
            'namd2':  NAMD2(self.config),
            'data':   Filewriter()
        }
        task_classes=inspect_classes('pestifer.tasks')
        self.tasks=[]
        prior_task=None
        for taskdict in self.config.tasks:
            assert len(taskdict)==1
            taskname=list(taskdict.keys())[0]
            specval=list(taskdict.values())[0]
            specs={} if not specval else specval.copy()
            class_name=[name for name,cls in task_classes.items() if cls.yaml_header==taskname][0]
            Cls=task_classes[class_name]
            this_task=Cls(specs,taskname,self.config,self.writers,prior_task)
            self.tasks.append(this_task)
            prior_task=this_task
        if len(self.tasks)==0 or not self.tasks[-1].taskname=='terminate':
            logger.debug('Adding default terminate task')
            self.tasks.append(TerminateTask({},'terminate',self.config,self.writers,prior_task))

        logger.info(f'Controller will execute {len(self.tasks)} task(s).')

    def do_tasks(self):
        for task in self.tasks:
            task.do()
