"""

.. module:: config
   :synopsis: Manages all configuration input and parsing
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import logging
logger=logging.getLogger(__name__)
import yaml
from .resourcemanager import ResourcesSetup, ResourcesApplyUserOptions, ResourcesGet

def replace(data,match,repl):
    """Recursive value search-and-replace; data is either list or dictionary; nesting is ok

    :param data: list or dict
    :type data: list or dict
    :param match: replacement string; expected to be found as $(val) in values
    :type match: str
    :param repl: replacement
    :type repl: *
    """
    match_str=r'$('+match+r')'
    if isinstance(data,(dict,list)):
        for k,v in (data.items() if isinstance(data,dict) else enumerate(data)):
            if v==match_str:
                data[k]=repl
            elif type(v)==str and match_str in v:
                data[k]=data[k].replace(match_str,repl)
            replace(v,match,repl)

class Config:
    def __init__(self):
        pass
    def setup(self,userconfigfilename):
        self.resman=ResourcesSetup()
        defaultconfigfilename=os.path.join(self.resman.ResourcePaths['config'],'defaults.yaml')
        """ Read the system default config file into the class defs """
        with open(defaultconfigfilename,'r') as f:
            self.defs=yaml.safe_load(f)
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
            self.defs['Segname_Second_Character']={'PROTEIN':'','ION':'I','WATER':'W','GLYCAN':'G','LIGAND':'L','OTHER':'O'}
            user_defs={}
            for v,r in self.defs.get('User_defaults',{}).items():
                user_defs[v]=r

        """ Read the user-specified config file """
        if userconfigfilename:
           with open(userconfigfilename,'r') as f:
                # read into instance not class
                user_defs.update(yaml.safe_load(f))
        ResourcesApplyUserOptions(user_defs)
        """ Read any mods or modfiles """
        for step in user_defs.get('BuildSteps',[]):
            if 'mods' in step: # either a dictionary or a file containing a dictionary
                val=step['mods']
                if type(val)==str and (val.endswith('yaml') or val.endswith('yml')):
                    with open(val,"r") as f:
                        step['mods']=yaml.safe_load(f)
        self.defs.update(user_defs)
        """ Perform any variable substitions at leaves """
        vars={'HOME':ResourcesGet('UserHome')}
        for v,r in vars.items():
            replace(self.defs,v,r)
        
    def __str__(self):
        return yaml.dump(self.data)        

TheConfig=Config()

def ConfigSetup(userconfigfilename):
    TheConfig.setup(userconfigfilename)
    return TheConfig

def ConfigGetParam(key,subkey=''):
    param=TheConfig.defs[key]
    if subkey and type(param)==dict:
        param=param[subkey]
    return param