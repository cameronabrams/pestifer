"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)
from .command import Command
import os
from .resources import ResourceManager

class Psfgen():
    def __init__(self,resman:ResourceManager,**options):
        if not resman:
            logger.fatal(f'Cannot use psfgen without a resource manager.')
        self.script_name=options.get('psfgen_script_name','my_mkpsf.tcl')
        self.script=['#### BEGIN HEADER ####']
        for t in options['StdCharmmTopo']:
            ft=os.path.join(resman.system_charmm_toppardir,t)
            self.script.append(f'topology {ft}')
        for lt in options['LocalCharmmTopo']:
            ft=os.path.join(resman.resource_paths['charmm'],lt)
            self.script.append(f'topology {ft}')

    '''under construction'''

    def run(self):
        c=Command()
        c.run()