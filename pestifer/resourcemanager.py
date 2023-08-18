"""

.. module:: resources
   :synopsis: Manages resources needed for pestifer, including those embedded in the Resources subpackage as well as system resources
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import os
import platform
import glob
import shutil
def is_tool(name):
    return shutil.which(name) is not None

logger=logging.getLogger(__name__)

from .util import replace
from pestifer import Resources # subpackage data

excludes=['__pycache__','__init__.py']

class ResourceManager:
    def __init__(self):
        self.root=os.path.dirname(Resources.__file__)
        ResourceFullPaths=glob.glob(self.Root+'/*')
        self.resource_paths={}
        for l in ResourceFullPaths:
            bn=os.path.basename(l)
            if not bn in excludes:
                self.resource_paths[bn]=l
        self.platform=platform.system()
        if self.platform=='Linux':
            self.user_home=os.environ['HOME']
        elif self.platform=='Windows':
            self.user_home=os.environ['HOMEPATH']
        elif self.platform=='Darwin':
            self.user_home=os.environ['HOME']
        else:
            raise(Exception,f'platform {self.platform} not recognized')
        
    def ApplyUserOptions(self,useroptions):
        replace(useroptions,'$(HOME)',self.user_home)

        self.user_charmmpath=useroptions.get('CHARMMPATH','')
        if not os.path.isdir(self.user_charmmpath):
            logger.warning(f'Your configuration indicates the charmm force-field files are located in {self.user_charmmpath}, but this directory is not found.')
        charmm_topparpath=os.path.join(self.user_charmmpath,'toppar')
        self.user_charmm_topparpath=None
        if os.path.exists(charmm_topparpath):
            self.system_charmm_topparpath=charmm_topparpath
        else:
            logger.warning(f'No directory "toppar/" detected in {self.user_charmmpath}')

        self.namd2=useroptions.get('NAMD2','')
        if not is_tool(self.namd2):
            logger.warning(f'{self.namd2}: Not found.')
        self.charmrun=useroptions.get('CHARMRUN','')
        if not is_tool(self.charmrun):
            logger.warning(f'{self.charmrun}: Not found.')        
        self.vmd=useroptions.get('VMD','')
        if not is_tool(self.vmd):
            logger.warning(f'{self.vmd}: Not found.')
        self.vmd_startup_script=useroptions.get('vmd_startup','pestifer-vmd.tcl')
        self.default_vmd_script=useroptions.get('default_vmd_script','vmd.tcl')
        self.default_psfgens_script=useroptions.get('default_psfge_script','psfgen.tcl')
        self.default_namd2_config=useroptions.get('default_namd2_config','namd2.tcl')

        return self

    def __str__(self):
        msg=f'Pestifer resources are located in "{self.root}"\n'
        for k,v in self.resource_paths.items():
            msg+=f'    {k:>20s}: {v}'
        msg+=f'External CHARMM force-field files at {self.user_charmm_topparpath}'
        msg+=f'charmrun executable as "{self.charmrun}"'
        msg+=f'namd2 executable as "{self.namd2}"'
        msg+=f'vmd executable as "{self.vmd}"'
        return msg
