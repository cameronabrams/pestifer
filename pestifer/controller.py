"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

from .config import ConfigSetup
from .steptask import Step
class Controller:
    def __init__(self,userconfigfilename):
        self.config=ConfigSetup(userconfigfilename)
        self.steps=[]
        self.register_mod_classes()
        self.register_modlist_classes()
        if 'steps' in self.config.defs:
            for step in self.config.defs['steps']:
                self.steps.append(Step(step).resolve_tasks())

    def do_steps(self,**kwargs):
        self.check()
        for step in self.steps:
            step.injest_molecules()
            step.do_tasks()

    def check(self):
        stepnames=list(set([x.name for x in self.steps]))
        assert len(stepnames)==len(self.steps),f'Please use unique step names in your "steps" section'
