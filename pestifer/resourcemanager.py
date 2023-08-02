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
from pestifer import Resources # subpackage data
# print(pkg_resources.files(Resources))

excludes=['__pycache__','__init__.py']

class ResourceManager:
    # def __init__(ResourceManager):
    # fullname=Resources.__file__
    ResourcesRoot=os.path.dirname(Resources.__file__)
    ResourceFullPaths=glob.glob(ResourcesRoot+'/*')
    ResourcePaths={}
    for l in ResourceFullPaths:
        bn=os.path.basename(l)
        if not bn in excludes:
            ResourcePaths[bn]=l
    System={}
    System['platform']=platform.system()
    if System['platform']=='Linux':
        System['user_home']=os.environ['HOME']
    elif System['platform']=='Windows':
        System['user_home']=os.environ['HOMEPATH']
    elif System['platform']=='Darwin':
        System['user_home']=os.environ['HOME']
    else:
        raise(Exception,f'platform {System["platform"]} not recognized')
    @classmethod
    def ApplyUserOptions(ResourceManager,useroptions):
        ResourceManager.system_charmmdir=useroptions.get('CHARMMDIR',None)
        # print(ResourceManager.system_charmmdir)
        if not os.path.isdir(ResourceManager.system_charmmdir):
            logger.warning(f'Your configuration indicates the charmm force-field files are located in {ResourceManager.system_charmmdir}, but this directory is not found.')
        charmm_toppardir=os.path.join(ResourceManager.system_charmmdir,'toppar')
        ResourceManager.system_charmm_toppardir=None
        if os.path.isdir(charmm_toppardir):
            ResourceManager.system_charmm_toppardir=charmm_toppardir
        else:
            logger.warning(f'No directory "toppar/" detected in {ResourceManager.system_charmmdir}')

        ResourceManager.namd2=useroptions.get('NAMD2',None)
        if not os.path.isfile(ResourceManager.namd2):
            logger.warning(f'{ResourceManager.namd2}: Not found.')
        ResourceManager.charmrun=useroptions.get('CHARMRUN',None)
        if not os.path.isfile(ResourceManager.charmrun):
            logger.warning(f'{ResourceManager.charmrun}: Not found.')        
        ResourceManager.vmd=useroptions.get('VMD',None)
        if not os.path.isfile(ResourceManager.vmd):
            logger.warning(f'{ResourceManager.vmd}: Not found.')

    def __str__():
        msg=f'Resources root is "{ResourceManager.ResourcesRoot}"\n'
        for k,v in ResourceManager.ResourcePaths.items():
            msg+=f'    {k:>20s}: {v}'
        msg+=f'System-wide CHARMM force-field files at {ResourceManager.system_charmm_toppardir}'
        msg+=f'CHARMRUN at {ResourceManager.charmrun}'
        msg+=f'NAMD2 at {ResourceManager.namd2}'
        msg+=f'VMD at {ResourceManager.vmd}'
        return msg
