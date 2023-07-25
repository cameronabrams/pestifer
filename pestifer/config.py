"""

.. module:: config
   :synopsis: Manages all configuration input and parsing
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import logging
logger=logging.getLogger(__name__)
import yaml
from .resources import ResourceManager

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
    def __init__(self,userconfigfilename,resman:ResourceManager):
        defaultconfigfilename=os.path.join(resman.resource_paths['config'],'defaults.yaml')
        """ Read the system default config file """
        with open(defaultconfigfilename,'r') as f:
            self.data=yaml.safe_load(f)
            ### generate reversible dictionaries
            self.data['CHARMM_to_PDB_Resnames']={v:k for k,v in self.data['PDB_to_CHARMM_Resnames'].items()}
            self.data['PDB_3char_to_1char_Resnames']={v:k for k,v in self.data['PDB_1char_to_3char_Resnames'].items()}
            self.data['Segtypes_by_Resnames']={'HOH':'WATER'}
            self.data['Segtypes_by_Resnames'].update({k:'ION' for k in self.data['PDB_Ions']})
            self.data['Segtypes_by_Resnames'].update({k:'GLYCAN' for k in self.data['PDB_Glycans']})
            self.data['Segtypes_by_Resnames'].update({self.data['PDB_to_CHARMM_Resnames'][k]:'GLYCAN' for k in self.data['PDB_Glycans']})
            self.data['Segtypes_by_Resnames'].update({k:'LIGAND' for k in self.data['PDB_Ligands']})
            self.data['Segtypes_by_Resnames'].update({self.data['PDB_to_CHARMM_Resnames'][k]:'LIGAND' for k in self.data['PDB_Ligands']})
            self.data['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in self.data['PDB_1char_to_3char_Resnames'].values()})
            self.data['Segtypes_by_Resnames'].update({k:'PROTEIN' for k in ['HIS','HSE','HSD']})
            self.data['Segname_Second_Character']={'PROTEIN':'','ION':'I','WATER':'W','GLYCAN':'G','LIGAND':'L','OTHER':'O'}
        """ Read the user-specified config file """
        with open(userconfigfilename,'r') as f:
            self.data.update(yaml.safe_load(f))
        """ Read any mods or modfiles """
        if 'BuildSteps' in self.data:
            for step in self.data['BuildSteps']:
                if 'mods' in step: # either a dictionary or a file containing a dictionary
                    val=step['mods']
                    if type(val)==str and (val.endswith('yaml') or val.endswith('yml')):
                        with open(val,"r") as f:
                            step['mods']=yaml.safe_load(f)
        """ Perform any variable substitions at leaves """
        vars={'HOME':resman.user_home}
        for v,r in vars.items():
            replace(self.data,v,r)
    def __str__(self):
        return yaml.dump(self.data)        
    