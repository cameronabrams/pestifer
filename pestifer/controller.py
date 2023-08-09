"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

from .config import ConfigSetup
from .molecule import Molecule

class Controller:
    def __init__(self,userconfigfilename):
        self.config=ConfigSetup(userconfigfilename)
    
    def build_molecules(self):
        self.molecules=[]
        if 'SourcePDB' in self.config.data:
            self.molecules.append(Molecule.from_pdb(pdb_code=self.config.data['SourcePDB']))

    # def report(self):
    #     print(str(self.config))
    #     print(str(self.resman))