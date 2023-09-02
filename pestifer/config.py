"""

.. module:: config
   :synopsis: Manages all configuration input and parsing, and access to installed resources
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import logging
logger=logging.getLogger(__name__)
import yaml
import glob
from collections import UserDict
# from .util import special_update, replace, is_tool
from pestifer import PestiferResources

class ResourceManager(UserDict):
    excludes=['__pycache__','_archive','bash']
    def __init__(self):
        data={}
        data['root']=os.path.dirname(PestiferResources.__file__)
        roottk=data['root'].split(os.sep)
        rootdepth=len(roottk)
        for root,dirs,files in os.walk(data['root']):
            tk=root.split(os.sep)
            absdepth=len(tk)
            if rootdepth<absdepth<rootdepth+3:
                rn=tk[rootdepth:absdepth]
                bn=tk[-1]
                if not bn in ResourceManager.excludes:
                    # logger.debug(f'{rn} {root}')
                    myset(data,rn,root)
                    
        # self._walk(data,[x for x in glob.glob(data['root']+'/*') if os.path.isdir(x)])
        super().__init__(data)

def myset(a_dict,keylist,val):
    if len(keylist)==1:
        a_dict[keylist[0]]=val
    else:
        key=keylist.pop(0)
        if not key in a_dict or type(a_dict[key])==type(val):
            a_dict[key]={}
        myset(a_dict[key],keylist,val)

class Config(UserDict):
    def __init__(self,userconfigfile=''):
        data={}
        data['Resources']=ResourceManager()
        logger.debug(f'Resources {data["Resources"]}')
        l=glob.glob(data['Resources']['config']+'/*.yaml')
        logger.debug(f'{l}')
        for fn in l:
            with open(fn,'r') as f:
                bn=os.path.basename(fn)
                bn=os.path.splitext(bn)[0]
                logger.debug(f'Pestifer uses {bn}: {fn}')
                data[bn]=yaml.safe_load(f)
        assert 'base' in data
        assert 'userhelp' in data
        if userconfigfile:
            if os.path.exists(userconfigfile):
                logger.debug(f'Pestifer uses {userconfigfile}: {fn}')
                with open(userconfigfile,'r') as f:
                    data['user']=yaml.safe_load(f)
                    self._user_defaults()
        super().__init__(data)
        self._set_shortcuts()

    def _user_defaults(self):
        # top-level directives, except for tasks
        for k,sp in self['userhelp'].items():
            if k!='tasks':
                if 'default' in sp.keys():
                    dv=sp['default']
                    if not k in self['user']:
                        self['user'][k]=dv


    def _set_shortcuts(self):
        self.tcl_proc_path=self['Resources']['tcl']['proc']
        self.tcl_script_path=self['Resources']['tcl']['scripts']
        self.vmd_startup_script=os.path.join(self.tcl_script_path,'pestifer-vmd.tcl')
        self.charmmff_toppar_path=self['Resources']['charmmff']['toppar']
        self.charmmff_custom_path=self['Resources']['charmmff']['custom']
        self.user_charmmff_toppar_path=''
        if hasattr(self,'user'):
            self.user_charmmff_toppar_path=os.path.join(self['user']['charmff'],'toppar')
        self.namd2_config_defaults=self['base']['namd2']
        for stn,stspec in self['base']['psfgen']['segtypes'].items():
            if stspec and 'rescodes' in stspec:
                self['base']['psfgen']['segtypes'][stn]['invrescodes']={v:k for k,v in stspec['rescodes'].items()}
        self.pdb_to_charmm_resnames={}
        for alias in self['base']['psfgen']['aliases']:
            tok=alias.split()
            if tok[0]=='residue':
                self.pdb_to_charmm_resnames[tok[1]]=tok[2]
    
    def segtype_by_resname(self,resname):
        sts=self['base']['psfgen']['segtypes']
        for s,rl in sts.items():
            if resname in rl['resnames']:
                return s
        return 'other'
    
    def charmmify_resname(self,pdbname):
        return self.pdb_to_charmm_resnames.get(pdbname,pdbname)
    
    def res_321(self,resname):
        return self['base']['psfgen']['segtypes']['protein']['rescodes'][resname]

    def res_123(self,rescode):
        return self['base']['psfgen']['segtypes']['protein']['invrescodes'][rescode]

#         self.charmm_custom_path=self.resman.resource_paths['charmm']['custom']
#         self.user_charmm_path=self.get('outside_charmm_path','UNSET')
#         self.user_charmm_toppar_path=os.path.join(self.user_charmm_path,'toppar')
# class Config(UserDict):
#     executables=['namd2','charmrun','vmd']
#     def __init__(self,*objs):
#         self.resman=ResourceManager()
#         data=self._readmainconfig()
#         if len(objs)==1:
#             userconfigfilename=objs[0]
#         else:
#             userconfigfilename=None
#         self._readuserconfig(userconfigfilename,data)
#         self._builddicts(data)
#         vars={'HOME':self.resman.user_home}
#         for v,r in vars.items():
#             replace(data,v,r)
#         super().__init__(data)
#         self.direct_attributes()

#     def direct_attributes(self):
#         self.has_executable={} 
#         for x in self.executables:
#             if not is_tool(self[x]):
#                 logger.warning(f'Executable "{x}" is not found.')
#                 self.has_executable[x]=False
#             else:
#                 self.has_executable[x]=True
#         self.vmd_startup_script=os.path.join(self.resman.resource_paths['tcl'],self.get('vmd_startup','pestifer-vmd.tcl'))
#         self.tcl_path=self.resman.resource_paths['tcl']
#         self.charmm_toppar_path=self.resman.resource_paths['charmm']['toppar']
#         self.charmm_custom_path=self.resman.resource_paths['charmm']['custom']
#         self.user_charmm_path=self.get('outside_charmm_path','UNSET')
#         self.user_charmm_toppar_path=os.path.join(self.user_charmm_path,'toppar')
#         if os.path.exists(self.user_charmm_path) and not os.path.exists(self.user_charmm_toppar_path):
#             logger.warning(f'You have specified the location where you keep your preferred CHARMM force-field/topology files, but no toppar/ directory found there.')

#     def _readmainconfig(self):
#         mainconfigfilename=os.path.join(self.resman.resource_paths['config'],'defaults.yaml')
#         assert os.path.exists(mainconfigfilename),f'Main configuration file {mainconfigfilename} not found.  Your installation of pestifer is corrupt.'
#         with open(mainconfigfilename,'r') as f:
#             data=yaml.safe_load(f)
#         self.defaults=data['User_defaults']
#         self.namd_params=data['Namd_defaults']
#         del data['User_defaults']
#         del data['Namd_defaults']
#         return data

#     def _readuserconfig(self,userconfigfilename,data):
#         user_dict={}
#         if userconfigfilename:
#            with open(userconfigfilename,'r') as f:
#                 user_dict=yaml.safe_load(f)
#         self.tasks=user_dict.get('tasks',{})
#         if 'tasks' in user_dict:
#             del user_dict['tasks']
#         namd_params=user_dict.get('Namd_params',{})
#         special_update(self.namd_params,namd_params)
#         if namd_params:
#             del user_dict['Namd_params']
#         tmp=self.defaults.copy()
#         tmp.update(user_dict)
#         special_update(data,tmp)

#     def segtype_resname_map(self):
#         return self['Segtypes_by_Resnames']
    
#     def segtype(self,resname):
#         return self['Segtypes_by_Resnames'].get(resname,None)

#     def _builddicts(self,data):
#         data['PDB_to_CHARMM_Resnames']={}
#         data['CHARMM_to_PDB_Resnames']={}
#         for a in data['Psfgen_Aliases']:
#             tok=a.split()
#             if tok[0]=='residue':
#                 data['PDB_to_CHARMM_Resnames'][tok[1]]=tok[2]
#                 data['CHARMM_to_PDB_Resnames'][tok[2]]=tok[1]
#         data['PDB_3char_to_1char_Resnames']={v:k for k,v in data['PDB_1char_to_3char_Resnames'].items()}
#         data['Segtypes_by_Resnames']={}
#         for segtype in data['Segtypes']:
#             resnames=data.get(segtype,[])
#             for r in resnames:
#                 data['Segtypes_by_Resnames'][r]=segtype
#                 if r in data['PDB_to_CHARMM_Resnames']:
#                     data['Segtypes_by_Resnames'][data['PDB_to_CHARMM_Resnames'][r]]=segtype
#         data['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in ['HIS','HSE','HSD']})

#     def __str__(self):
#         return yaml.dump(self)        
