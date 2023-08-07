"""

.. module:: resources
   :synopsis: Manages resources needed for pestifer, including those embedded in the Resources subpackage as well as system resources
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import os
import platform
import glob
logger=logging.getLogger(__name__)

from pestifer import Resources # subpackage data

excludes=['__pycache__','__init__.py']

class ResourceManager:
    def __init__(self):
        self.Root=os.path.dirname(Resources.__file__)

    def setup(self):
        ResourceFullPaths=glob.glob(self.Root+'/*')
        self.ResourcePaths={}
        for l in ResourceFullPaths:
            bn=os.path.basename(l)
            if not bn in excludes:
                self.ResourcePaths[bn]=l
        self.Platform=platform.system()
        if self.Platform=='Linux':
            self.UserHome=os.environ['HOME']
        elif self.Platform=='Windows':
            self.UserHome=os.environ['HOMEPATH']
        elif self.Platform=='Darwin':
            self.UserHome=os.environ['HOME']
        else:
            raise(Exception,f'platform {self.Platform} not recognized')

    def ApplyUserOptions(self,useroptions):
        self.system_charmmdir=useroptions.get('CHARMMDIR','')
        # print(self.system_charmmdir)
        if not os.path.isdir(self.system_charmmdir):
            logger.warning(f'Your configuration indicates the charmm force-field files are located in {self.system_charmmdir}, but this directory is not found.')
        charmm_toppardir=os.path.join(self.system_charmmdir,'toppar')
        self.system_charmm_toppardir=None
        if os.path.exists(charmm_toppardir):
            self.system_charmm_toppardir=charmm_toppardir
        else:
            logger.warning(f'No directory "toppar/" detected in {self.system_charmmdir}')

        self.namd2=useroptions.get('NAMD2','')
        if not os.path.isfile(self.namd2):
            logger.warning(f'{self.namd2}: Not found.')
        self.charmrun=useroptions.get('CHARMRUN','')
        if not os.path.isfile(self.charmrun):
            logger.warning(f'{self.charmrun}: Not found.')        
        self.vmd=useroptions.get('VMD','')
        if not os.path.isfile(self.vmd):
            logger.warning(f'{self.vmd}: Not found.')

    def __str__(self):
        msg=f'Resources root is "{self.Root}"\n'
        for k,v in self.ResourcePaths.items():
            msg+=f'    {k:>20s}: {v}'
        msg+=f'System-wide CHARMM force-field files at {self.system_charmm_toppardir}'
        msg+=f'CHARMRUN at {self.charmrun}'
        msg+=f'NAMD2 at {self.namd2}'
        msg+=f'VMD at {self.vmd}'
        return msg

TheResourceManager=ResourceManager()

def ResourcesSetup():
    TheResourceManager.setup()
    return TheResourceManager

def ResourcesApplyUserOptions(useroptions):
    TheResourceManager.ApplyUserOptions(useroptions)
    return TheResourceManager

def ResourcesGet(attr):
    return TheResourceManager.__dict__.get(attr,None)

def ResourcesGetPath(key):
    return TheResourceManager.ResourcePaths.get(key,None)

import yaml
def ResourcesInfo():
    msg='# YAML dump of ResourceManager\n'
    msg+=yaml.dump(TheResourceManager.__dict__)
    return msg