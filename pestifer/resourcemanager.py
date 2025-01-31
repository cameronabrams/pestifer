# Author: Cameron F. Abrams, <cfa22@drexel.edu>
import glob
import os
from . import resources
import logging
logger=logging.getLogger(__name__)
class ResourceManager:
    base_resources=['charmmff','examples','tcl','ycleptic']
    ignored_resources=['__pycache__','_archive','bash']
    def __init__(self):
        self.pdb_depot_path=None
        self.resources_path=os.path.dirname(resources.__file__)
        self.resource_dirs=[x for x in glob.glob(os.path.join(self.resources_path,'*')) if os.path.isdir(x) and not os.path.basename(x) in ResourceManager.ignored_resources]
        # self.resource_dirs=[x for x in glob.glob(os.path.join(self.resources_path,'*')) if not os.path.basename(x) in ResourceManager.ignored_resources]
        # assert all([x in self.resource_dirs for x in ResourceManager.base_resources]),'Error: this installation of pestifer is missing Resources.'
        assert all([x in [os.path.basename(_) for _ in self.resource_dirs] for x in ResourceManager.base_resources]),f'some resources seem to be missing'
        self.ycleptic_configdir=os.path.join(self.resources_path,'ycleptic')
        ycleptic_files=glob.glob(os.path.join(self.ycleptic_configdir,'*'))
        assert len(ycleptic_files)==1,f'Too many config files in {self.ycleptic_configdir}: {ycleptic_files}'
        self.ycleptic_config=ycleptic_files[0]
        self.resource_path={}
        for r in ResourceManager.base_resources:
            self.resource_path[r]=os.path.join(self.resources_path,r)
    
    def __str__(self):
        retstr='Pestifer Resources are found in these package subdirectories:'
        for r,p in self.resource_path.items():
            retstr+=f'{p}\n'
        return retstr

    def get_ycleptic_config(self):
        return self.ycleptic_config
    
    def get_resource_path(self,r):
        return self.resource_path.get(r,None)
    
    def get_examples_as_list(self,fullpaths=False):
        epath=self.resource_path['examples']
        fullnames=glob.glob(os.path.join(epath,'*'))
        fullnames.sort()
        if fullpaths:
            return fullnames
        basenames=[os.path.basename(x) for x in fullnames]
        basenames.sort()
        return basenames
    
    def get_example_yaml_by_index(self,index):
        epath=self.resource_path['examples']
        basenames=self.get_examples_as_list()
        logger.debug(epath)
        logger.debug(epath)
        sindex=f'{index:02d}'
        for b in basenames:
            if b.startswith(sindex):
                break
        else:
            return None
        return os.path.join(epath,b)

    def get_charmmff_toppardir(self):
        return os.path.join(self.resource_path['charmmff'],'toppar')
    
    def get_charmmff_customdir(self):
        return os.path.join(self.resource_path['charmmff'],'custom')

    def get_charmmff_pdbdir(self):
        ### pestifer's installed pdb depot
        return os.path.join(self.resource_path['charmmff'],'pdb')
    
    def get_charmmff_pdbdepot(self):
        ### user-specified pdb depot
        return self.pdb_depot_path

    def get_charmmff_pdb_streams_as_list(self,fullpaths=False):
        pdbdir=self.get_charmmff_pdbdir()
        streams_full=[x for x in glob.glob(os.path.join(pdbdir,'*')) if os.path.isdir(x)]
        if self.pdb_depot_path:
            user_streams_full=[x for x in glob.glob(os.path.join(self.pdb_depot_path,'*')) if os.path.isdir(x)]
            streams_full.extend(user_streams_full)
        streams=[os.path.basename(x) for x in streams_full]
        return streams
    
    def get_charmmff_pdbs_as_dict_by_stream(self):
        pdbdir=self.get_charmmff_pdbdir()
        streams=self.get_charmmff_pdb_streams_as_list()
        pdbdirs_dict={}
        for stream in streams:
            streamdir=os.path.join(pdbdir,stream)
            pdbdirs_full=[x for x in glob.glob(os.path.join(streamdir,'*')) if os.path.isdir(x)]
            raw_pdbs=[os.path.splitext(os.path.basename(x))[0] for x in glob.glob(os.path.join(streamdir,'*')) if x.endswith('.pdb')]
            raw_pdbs.sort()
            pdbdirs=[os.path.basename(x) for x in pdbdirs_full]
            pdbdirs.sort()
            pdbdirs_dict[stream]=pdbdirs if len(pdbdirs)>0 else raw_pdbs
        return pdbdirs_dict
    
    def get_charmmff_pdb_path(self,name):
        D=self.get_charmmff_pdbs_as_dict_by_stream()
        logger.debug(f'{D}')
        for s,L in D.items():
            if name in L:
                return os.path.join(self.get_charmmff_pdbdir(),s,name)
        return None        

    def get_tcldir(self):
        return self.resource_path['tcl']
    
    def get_tcl_pkgdir(self):
        return os.path.join(self.resource_path['tcl'],'pkg')
    
    def get_tcl_scriptsdir(self):
        return os.path.join(self.resource_path['tcl'],'scripts')
    
    def set_pdb_depot(self,depot_path):
        self.pdb_depot_path=depot_path