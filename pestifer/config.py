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

class Config(dict):
    executables=['namd2','charmrun','vmd']
    def __init__(self,userconfigfilename):
        self.resman=ResourceManager()
        self._readmainconfig()
        self._readuserconfig(userconfigfilename)
        # append/update any container-like values, rather than overwriting
        self.defs=special_update(self.defs,self.user_defaults)
        self._builddicts()

        vars={'HOME':self.resman.user_home}
        for v,r in vars.items():
            replace(self.defs,v,r)
        self.has_executable={} 
        for x in self.executables:
            tool=self.defs[x]
            if not is_tool(tool):
                logger.warning(f'Executable "{x}" is not found.')
                self.has_executable[x]=False
            else:
                self.has_executable[x]=True

        self.vmd_startup_script=os.path.join(self.resman.resource_paths['tcl'],self.defs.get('vmd_startup','pestifer-vmd.tcl'))
        self.tcl_path=self.resman.resource_paths['tcl']
        self.namd_template_path=self.resman.resource_paths['templates']
        self.charmm_toppar_path=self.resman.resource_paths['charmm']
        self.user_charmm_path=self.defs['charmm_path']
        self.user_charmm_toppar_path=os.path.join(self.user_charmm_path,'toppar')
        self.has_charmm=True
        if not os.path.exists(self.user_charmm_toppar_path):
            logger.warning(f'You have not specified the location where you keep your CHARMM force-field/topology files.')
            logger.warning(f'Also, the default location where I expect to find them, {self.user_defaults["charmm_path"]}, seems not to exist.')
            self.has_charmm=False

    def _readmainconfig(self):
        mainconfigfilename=os.path.join(self.resman.resource_paths['config'],'defaults.yaml')
        assert os.path.exists(mainconfigfilename),f'Main configuration file {mainconfigfilename} not found.  Your installation of pestifer is corrupt.'
        with open(mainconfigfilename,'r') as f:
            self.defs=yaml.safe_load(f)
        self.user_defaults=self.defs.get('User_defaults',{})

    def _readuserconfig(self,userconfigfilename):
        if userconfigfilename:
           with open(userconfigfilename,'r') as f:
                self.user_dict=yaml.safe_load(f)
                tmp=self.user_defaults.copy()
                tmp.update(self.user_dict)
                self.user_dict=tmp

    def __getitem__(self,key):
        return self.defs[key]
    
    def get(self,key,default):
        if key in self.defs:
            return self.defs[key]
        else:
            return default

    def segtype_resname_map(self):
        return self.defs['Segtypes_by_Resnames']
    
    def segtype(self,resname):
        return self.defs['Segtypes_by_Resnames'].get(resname,None)

    def _builddicts(self):
        self.defs['CHARMM_to_PDB_Resnames']={v:k for k,v in self.defs['PDB_to_CHARMM_Resnames'].items()}
        self.defs['PDB_3char_to_1char_Resnames']={v:k for k,v in self.defs['PDB_1char_to_3char_Resnames'].items()}
        self.defs['Segtypes_by_Resnames']={'HOH':'WATER'}
        self.defs['Segtypes_by_Resnames'].update({k:'ION' for k in self.defs['PDB_Ions']})
        self.defs['Segtypes_by_Resnames'].update({k:'GLYCAN' for k in self.defs['PDB_Glycans']})
        self.defs['Segtypes_by_Resnames'].update({self.defs['PDB_to_CHARMM_Resnames'][k]:'GLYCAN' for k in self.defs['PDB_Glycans']})
        self.defs['Segtypes_by_Resnames'].update({k:'LIGAND' for k in self.defs['PDB_Ligands']})
        self.defs['Segtypes_by_Resnames'].update({self.defs['PDB_to_CHARMM_Resnames'][k]:'LIGAND' for k in self.defs['PDB_Ligands']})
        self.defs['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in self.defs['PDB_1char_to_3char_Resnames'].values()})
        self.defs['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in ['HIS','HSE','HSD']})


    def __str__(self):
        return yaml.dump(self)        
