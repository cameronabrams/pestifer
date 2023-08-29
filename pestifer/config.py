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
from pestifer import Resources # subpackage data

class ResourceManager:
    excludes=['__pycache__','__init__.py']
    def __init__(self):
        self.root=os.path.dirname(Resources.__file__)
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
        self.charmm_toppar_path=self.resman.resource_paths['charmm']
        self.user_charmm_path=self['charmm_path']
        self.user_charmm_toppar_path=os.path.join(self.user_charmm_path,'toppar')
        self.has_charmm=True
        if not os.path.exists(self.user_charmm_toppar_path):
            logger.warning(f'You have not specified the location where you keep your CHARMM force-field/topology files.')
            logger.warning(f'Also, the default location where I expect to find them, {self.defaults["charmm_path"]}, seems not to exist.')
            self.has_charmm=False

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
        data['CHARMM_to_PDB_Resnames']={v:k for k,v in data['PDB_to_CHARMM_Resnames'].items()}
        data['PDB_3char_to_1char_Resnames']={v:k for k,v in data['PDB_1char_to_3char_Resnames'].items()}
        data['Segtypes_by_Resnames']={'HOH':'WATER'}
        data['Segtypes_by_Resnames'].update({k:'ION' for k in data['PDB_Ions']})
        data['Segtypes_by_Resnames'].update({k:'GLYCAN' for k in data['PDB_Glycans']})
        data['Segtypes_by_Resnames'].update({data['PDB_to_CHARMM_Resnames'][k]:'GLYCAN' for k in data['PDB_Glycans']})
        data['Segtypes_by_Resnames'].update({k:'LIGAND' for k in data['PDB_Ligands']})
        data['Segtypes_by_Resnames'].update({data['PDB_to_CHARMM_Resnames'][k]:'LIGAND' for k in data['PDB_Ligands']})
        data['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in data['PDB_1char_to_3char_Resnames'].values()})
        data['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in ['HIS','HSE','HSD']})

    def __str__(self):
        return yaml.dump(self)        
