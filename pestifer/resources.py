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

def replace(data, match, repl):
    if isinstance(data, (dict, list)):
        for k, v in (data.items() if isinstance(data, dict) else enumerate(data)):
            if v == match:
                data[k] = repl
            replace(v, match, repl)

class ResourceManager:
    def __init__(self,**kwargs):
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

    def ApplyUserOptions(self,**options):
        vars={'HOME':self.user_home}
        for v,r in vars.items():
            replace(options,r'$('+v+r')',r)
        
        self.system_charmmdir=options.get('CHARMDIR',None)
        if not os.path.isdir(self.system_charmmdir):
            logger.fatal('Your configuration indicates the charmm force-field files are located in {self.system_charmmdir}, but this directory is not found.')
        charmm_toppardir=os.path.join(self.system_charmmdir,'toppar')
        self.system_charmm_toppardir=None
        if os.path.isdir(charmm_toppardir):
            self.system_charmm_toppardir=charmm_toppardir
        else:
            logger.fatal(f'No directory "toppar" detected in {self.system_charmmdir}')

        namd2_path=options.get('NAMD2',None)
        if os.path.isfile(namd2_path):
            self.namd2=namd2_path
        else:
            logger.warning(f'No valid location for namd2 provided.')
        charmrun_path=options.get('CHARMRUN',None)
        if os.path.isfile(charmrun_path):
            self.charmrun=charmrun_path
        else:
            logger.warning(f'No valid location for charmrun provided.')
        
        vmd_path=options.get('VMD',None)
        if os.path.isfile(vmd_path):
            self.vmd_location=vmd_path
        else:
            logger.warning(f'No valid location for vmd provided.')
        
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
