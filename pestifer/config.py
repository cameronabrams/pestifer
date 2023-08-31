"""

.. module:: config
   :synopsis: Manages all configuration input and parsing, and access to installed resources
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import logging
logger=logging.getLogger(__name__)
import yaml
import platform
import glob
from collections import UserDict
from .util import special_update, replace, is_tool
import PestiferResources

class ResourceManager:
    excludes=['__pycache__','__init__.py']
    def __init__(self):
        self.root=os.path.dirname(PestiferResources.__file__)
        ResourceFullPaths=glob.glob(self.root+'/*')
        self.resource_paths={}
        for l in ResourceFullPaths:
            bn=os.path.basename(l)
            if not bn in ResourceManager.excludes:
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

        CharmmPaths=glob.glob(self.resource_paths['charmm']+'/*')
        self.resource_paths['charmm']={}
        for c in CharmmPaths:
            bn=os.path.basename(c)
            self.resource_paths['charmm'][bn]=c

    def __str__(self):
        msg=f'Pestifer resources are located in "{self.root}"\n'
        for k,v in self.resource_paths.items():
            msg+=f'    {k:>20s}: {v}'
        return msg

class Config(UserDict):
    executables=['namd2','charmrun','vmd']
    def __init__(self,*objs):
        self.resman=ResourceManager()
        data=self._readmainconfig()
        if len(objs)==1:
            userconfigfilename=objs[0]
        else:
            userconfigfilename=None
        self._readuserconfig(userconfigfilename,data)
        self._builddicts(data)
        vars={'HOME':self.resman.user_home}
        for v,r in vars.items():
            replace(data,v,r)
        super().__init__(data)
        self.direct_attributes()

    def direct_attributes(self):
        self.has_executable={} 
        for x in self.executables:
            if not is_tool(self[x]):
                logger.warning(f'Executable "{x}" is not found.')
                self.has_executable[x]=False
            else:
                self.has_executable[x]=True
        self.vmd_startup_script=os.path.join(self.resman.resource_paths['tcl'],self.get('vmd_startup','pestifer-vmd.tcl'))
        self.tcl_path=self.resman.resource_paths['tcl']
        self.charmm_toppar_path=self.resman.resource_paths['charmm']['toppar']
        self.charmm_custom_path=self.resman.resource_paths['charmm']['custom']
        self.user_charmm_path=self.get('outside_charmm_path','UNSET')
        self.user_charmm_toppar_path=os.path.join(self.user_charmm_path,'toppar')
        if os.path.exists(self.user_charmm_path) and not os.path.exists(self.user_charmm_toppar_path):
            logger.warning(f'You have specified the location where you keep your preferred CHARMM force-field/topology files, but no toppar/ directory found there.')

    def _readmainconfig(self):
        mainconfigfilename=os.path.join(self.resman.resource_paths['config'],'defaults.yaml')
        assert os.path.exists(mainconfigfilename),f'Main configuration file {mainconfigfilename} not found.  Your installation of pestifer is corrupt.'
        with open(mainconfigfilename,'r') as f:
            data=yaml.safe_load(f)
        self.defaults=data['User_defaults']
        self.namd_params=data['Namd_defaults']
        del data['User_defaults']
        del data['Namd_defaults']
        return data

    def _readuserconfig(self,userconfigfilename,data):
        user_dict={}
        if userconfigfilename:
           with open(userconfigfilename,'r') as f:
                user_dict=yaml.safe_load(f)
        self.tasks=user_dict.get('tasks',{})
        if 'tasks' in user_dict:
            del user_dict['tasks']
        namd_params=user_dict.get('Namd_params',{})
        special_update(self.namd_params,namd_params)
        if namd_params:
            del user_dict['Namd_params']
        tmp=self.defaults.copy()
        tmp.update(user_dict)
        special_update(data,tmp)

    def segtype_resname_map(self):
        return self['Segtypes_by_Resnames']
    
    def segtype(self,resname):
        return self['Segtypes_by_Resnames'].get(resname,None)

    def _builddicts(self,data):
        data['PDB_to_CHARMM_Resnames']={}
        data['CHARMM_to_PDB_Resnames']={}
        for a in data['Psfgen_Aliases']:
            tok=a.split()
            if tok[0]=='residue':
                data['PDB_to_CHARMM_Resnames'][tok[1]]=tok[2]
                data['CHARMM_to_PDB_Resnames'][tok[2]]=tok[1]
        data['PDB_3char_to_1char_Resnames']={v:k for k,v in data['PDB_1char_to_3char_Resnames'].items()}
        data['Segtypes_by_Resnames']={}
        for segtype in data['Segtypes']:
            resnames=data.get(segtype,[])
            for r in resnames:
                data['Segtypes_by_Resnames'][r]=segtype
                if r in data['PDB_to_CHARMM_Resnames']:
                    data['Segtypes_by_Resnames'][data['PDB_to_CHARMM_Resnames'][r]]=segtype
        data['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in ['HIS','HSE','HSD']})

    def __str__(self):
        return yaml.dump(self)        
