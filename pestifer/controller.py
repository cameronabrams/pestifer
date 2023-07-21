"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

from .resources import ResourceManager
from .config import Config

class Controller:
    def __init__(self,userconfigfilename):
        self.resman=ResourceManager()
        self.config=Config(userconfigfilename,self.resman)
        self.resman.ApplyUserOptions(self.config.data)
    
    def report(self):
        print(str(self.config))
        print(str(self.resman))