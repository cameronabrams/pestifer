"""

.. module:: resources
   :synopsis: Manages resources needed for pestifer, including those embedded in the Resources subpackage as well as system resources
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import os
import platform
import glob
from pathlib import Path
# import importlib.resources as pkg_resources
logger=logging.getLogger(__name__)
from pestifer import Resources
# print(pkg_resources.files(Resources))

excludes=['__pycache__','__init__.py']

class ResourceManager:
    def __init__(self):
        fullname=Resources.__file__
        self.ResourcesPath=os.path.dirname(fullname)
        l_elements=glob.glob(self.ResourcesPath+'/*')
        self.resource_paths={}
        for l in l_elements:
            bn=os.path.basename(l)
            if not bn in excludes:
                self.resource_paths[bn]=l
        list_dirs=os.walk(self.ResourcesPath)
        self.dirs=[]
        self.files=[]
        for root,dirs,files in list_dirs: 
            for d in dirs:
                pn=Path(os.path.join(root,d))
                excl=any([x in excludes for x in pn.parts])
                if not excl:
                    self.dirs.append(pn)
            for f in files:
                pn=Path(os.path.join(root,f))
                excl=any([x in excludes for x in pn.parts])
                if not excl:
                    self.files.append(pn)
        
        plat=platform.system()
        if plat=='Linux':
            self.user_home=os.environ['HOME']
        elif plat=='Windows':
            self.user_home=os.environ['HOMEPATH']
        elif plat=='Darwin':
            self.user_home=os.environ['HOME']
        else:
            raise(Exception,f'platform {plat} not recognized')

    def ApplyUserOptions(self,useroptions):
        self.system_charmmdir=useroptions.get('CHARMMDIR',None)
        # print(self.system_charmmdir)
        if not os.path.isdir(self.system_charmmdir):
            logger.warning(f'Your configuration indicates the charmm force-field files are located in {self.system_charmmdir}, but this directory is not found.')
        charmm_toppardir=os.path.join(self.system_charmmdir,'toppar')
        self.system_charmm_toppardir=None
        if os.path.isdir(charmm_toppardir):
            self.system_charmm_toppardir=charmm_toppardir
        else:
            logger.warning(f'No directory "toppar/" detected in {self.system_charmmdir}')

        self.namd2=useroptions.get('NAMD2',None)
        if not os.path.isfile(self.namd2):
            logger.warning(f'{self.namd2}: Not found.')
        self.charmrun=useroptions.get('CHARMRUN',None)
        if not os.path.isfile(self.charmrun):
            logger.warning(f'{self.charmrun}: Not found.')        
        self.vmd=useroptions.get('VMD',None)
        if not os.path.isfile(self.vmd):
            logger.warning(f'{self.vmd}: Not found.')
        
    def __str__(self):
        msg=f'Resources path is "{self.ResourcesPath}"\n'
        msg+=f'Resources directories:\n'
        for d in self.dirs:
            msg+=f'DIR  {d}\n'
        msg+=f'Resources files:\n'
        for f in self.files:
            msg+=f'FILE {f}\n'
        msg+=f'System-wide CHARMM force-field files at {self.system_charmm_toppardir}'
        # msg+=f'CHARMM force-field files at {self.charmm_location}'
        return msg
