"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

from .config import Config
from .scriptwriters import Filewriter,Psfgen,VMD,NAMD2
from .tasks import Task

class Controller:
    def __init__(self,userconfigfilename):
        self.config=Config(userconfigfilename)
        self.writers={
            'psfgen': Psfgen(self.config),
            'vmd':    VMD(self.config),
            'namd2':  NAMD2(self.config),
            'data':   Filewriter()
        }
        self.tasks=[]
        prior_task=None
        for taskdict in self.config.tasks:
            this_task=Task(taskdict,self.config,self.writers,prior_task)
            self.tasks.append(this_task)
            prior_task=this_task

        logger.debug(f'Controller will execute {len(self.tasks)} task(s).')

    def do_tasks(self):
        for task in self.tasks:
            task.do()
