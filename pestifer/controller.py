"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

from .config import ConfigSetup
from .psfgen import Psfgen,VMDScript
from .steptask import Step
class Controller:
    def __init__(self,userconfigfilename):
        self.config=ConfigSetup(userconfigfilename)
        self.psfgen=Psfgen(self.config.resman) # handler for psfgen scripts
        self.vmdtcl=VMDScript(self.config.resman) # handler for general (non-psfgen) vmd scripts
        self.steps=[]
        if 'steps' in self.config.defs:
            for stepdict in self.config.defs['steps']:
                self.steps.append(Step(stepdict).resolve_tasks(self.psfgen,self.vmdtcl))
        logger.debug(f'Controller will execute {len(self.steps)} step(s).')

    def do_steps(self,**kwargs):
        self.check()
        for step in self.steps:
            step.injest_molecules()
            step.do_tasks()

    def check(self):
        stepnames=list(set([x.name for x in self.steps]))
        assert len(stepnames)==len(self.steps),f'Please use unique step names in your "steps" section'
