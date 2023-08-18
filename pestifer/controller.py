"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

from .config import Config
from .scriptwriters import Psfgen,VMD,NAMD2
from .steptask import Step

class Controller:
    def __init__(self,userconfigfilename):
        self.config=Config(userconfigfilename)
        self.scriptwriters={
            'psfgen': Psfgen(self.config),
            'vmd':    VMD(self.config),
            'namd2':  NAMD2(self.config)
        }
        self.steps=[]
        if 'steps' in self.config.defs:
            for stepdict in self.config.defs['steps']:
                self.steps.append(Step(stepdict).resolve_tasks(self.config,self.scriptwriters))
        logger.debug(f'Controller will execute {len(self.steps)} step(s).')

    def do_steps(self,**kwargs):
        self._check()
        for step in self.steps:
            step.injest_molecules()
            step.do_tasks()

    def _check(self):
        stepnames=list(set([x.name for x in self.steps]))
        assert len(stepnames)==len(self.steps),f'Please use unique step names in your "steps" section'
