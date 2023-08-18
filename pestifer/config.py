"""

.. module:: config
   :synopsis: Manages all configuration input and parsing
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import logging
logger=logging.getLogger(__name__)
import yaml
from .util import special_update, replace

class Config(dict):
    def __init__(self,resource_manager,userconfigfilename):
        mainconfigfilename=os.path.join(resource_manager.resource_paths['config'],'defaults.yaml')
        assert os.path.exists(mainconfigfilename),f'Main configuration file {mainconfigfilename} not found.  Your installation of pestifer is corrupt.'

        """ Read the system default config file into the class defs """
        with open(mainconfigfilename,'r') as f:
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
            # self.defs['Segname_Second_Character']={'PROTEIN':'','ION':'I','WATER':'W','GLYCAN':'G','LIGAND':'L','OTHER':'O'}
            user_defs={}
            for v,r in self.defs.get('User_defaults',{}).items():
                user_defs[v]=r

        """ Read the user-specified config file """
        if userconfigfilename:
           with open(userconfigfilename,'r') as f:
                user_dict=yaml.safe_load(f)
                # do not allow 'User_defaults' key in a user-supplied yaml file
                if 'User_defaults' in user_dict:
                    logger.warning(f'Key "User_defaults" detected in user-supplied config file {userconfigfilename} is ignored.')
                    del user_dict['User_defaults']
                user_defs.update(user_dict)
        """ special_update appends to any list or dict values in the first arg """
        self.defs=special_update(self.defs,user_defs)
        """ Perform any variable substitions at leaves """
        vars={'HOME':resource_manager.user_home}
        for v,r in vars.items():
            replace(self.defs,v,r)

    def __getitem__(self,key):
        return self.defs[key]
    
    def __str__(self):
        return yaml.dump(self.data)        
