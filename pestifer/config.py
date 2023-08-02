"""

.. module:: config
   :synopsis: Manages all configuration input and parsing
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import logging
logger=logging.getLogger(__name__)
import yaml
from .resourcemanager import ResourceManager

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
    defs={}
    def __init__(self,userconfigfilename):
        self.user_data={}
        defaultconfigfilename=os.path.join(ResourceManager.ResourcePaths['config'],'defaults.yaml')
        """ Read the system default config file into the class defs """
        with open(defaultconfigfilename,'r') as f:
            Config.defs=yaml.safe_load(f)
            Config.defs['CHARMM_to_PDB_Resnames']={v:k for k,v in Config.defs['PDB_to_CHARMM_Resnames'].items()}
            Config.defs['PDB_3char_to_1char_Resnames']={v:k for k,v in Config.defs['PDB_1char_to_3char_Resnames'].items()}
            Config.defs['Segtypes_by_Resnames']={'HOH':'WATER'}
            Config.defs['Segtypes_by_Resnames'].update({k:'ION' for k in Config.defs['PDB_Ions']})
            Config.defs['Segtypes_by_Resnames'].update({k:'GLYCAN' for k in Config.defs['PDB_Glycans']})
            Config.defs['Segtypes_by_Resnames'].update({Config.defs['PDB_to_CHARMM_Resnames'][k]:'GLYCAN' for k in Config.defs['PDB_Glycans']})
            Config.defs['Segtypes_by_Resnames'].update({k:'LIGAND' for k in Config.defs['PDB_Ligands']})
            Config.defs['Segtypes_by_Resnames'].update({Config.defs['PDB_to_CHARMM_Resnames'][k]:'LIGAND' for k in Config.defs['PDB_Ligands']})
            Config.defs['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in Config.defs['PDB_1char_to_3char_Resnames'].values()})
            Config.defs['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in ['HIS','HSE','HSD']})
            Config.defs['Segname_Second_Character']={'PROTEIN':'','ION':'I','WATER':'W','GLYCAN':'G','LIGAND':'L','OTHER':'O'}
            for v,r in Config.defs.get('User_defaults',{}).items():
                self.user_data[v]=r

        """ Read the user-specified config file """
        if userconfigfilename:
           with open(userconfigfilename,'r') as f:
                # read into instance not class
                self.user_data.update(yaml.safe_load(f))
        ResourceManager.ApplyUserOptions(self.user_data)
        """ Read any mods or modfiles """
        for step in self.user_data.get('BuildSteps',[]):
            if 'mods' in step: # either a dictionary or a file containing a dictionary
                val=step['mods']
                if type(val)==str and (val.endswith('yaml') or val.endswith('yml')):
                    with open(val,"r") as f:
                        step['mods']=yaml.safe_load(f)
        """ Perform any variable substitions at leaves """
        vars={'HOME':ResourceManager.System['user_home']}
        for v,r in vars.items():
            replace(self.user_data,v,r)
            replace(Config.defs,v,r)

    # def __str__(self):
    #     return yaml.dump(self.data)        
    