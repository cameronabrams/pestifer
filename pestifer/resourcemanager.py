# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" Defines the ResourceManager class for managing access to pestifer's built-in resources 
"""
import glob
import os
from . import resources
from .pdbcollection import PDBCollection

import logging
logger=logging.getLogger(__name__)
class ResourceManager:
    base_resources=['charmmff','examples','tcl','ycleptic']
    ignored_resources=['__pycache__','_archive','bash']
    def __init__(self):
        self.resources_path=os.path.dirname(resources.__file__)
        self.resource_dirs=[x for x in glob.glob(os.path.join(self.resources_path,'*')) if os.path.isdir(x) and not os.path.basename(x) in ResourceManager.ignored_resources]
        assert all([x in [os.path.basename(_) for _ in self.resource_dirs] for x in ResourceManager.base_resources]),f'some resources seem to be missing'
        self.ycleptic_configdir=os.path.join(self.resources_path,'ycleptic')
        ycleptic_files=glob.glob(os.path.join(self.ycleptic_configdir,'*'))
        assert len(ycleptic_files)==1,f'Too many config files in {self.ycleptic_configdir}: {ycleptic_files}'
        self.ycleptic_config=ycleptic_files[0]
        self.resource_path={}
        for r in ResourceManager.base_resources:
            self.resource_path[r]=os.path.join(self.resources_path,r)
        self.pdb_collection=PDBCollection(os.path.join(self.resource_path['charmmff'],'pdb'))

    def __str__(self):
        cp=os.path.commonpath(list(self.resource_path.values()))
        retstr=f'Pestifer resources are found under\n    {cp}\n'
        for r,p in self.resource_path.items():
            retstr+=f'        {p.replace(cp+os.sep,"")+os.sep}\n'
        return retstr

    def show(self,out_stream=print,components={}):
        for c,spec in components.items():
            if not c in self.base_resources:
                logger.warning(f'{c} is not a base resource; expected one of {", ".join(self.base_resources)}')
            path=self.get_resource_path(c)
            if c=='examples':
                exb=self.get_examples_as_list()
                out_stream(f'\nExamples:')
                for e in exb:
                    out_stream(f'    {e}')
            elif c=='charmmff':
                if 'toppar' in spec:
                    path=self.get_charmmff_toppardir()
                    with open(os.path.join(path,'00PESTIFER-README.txt'),'r') as f:
                        msg=f.read()
                    out_stream(msg)
                if 'pdb' in spec:
                    self.pdb_collection.show(out_stream)
                if 'custom' in spec:
                    path=self.get_charmmff_customdir()
                    with open(os.path.join(path,'00PESTIFER-README.txt'),'r') as f:
                        msg=f.read()
                    out_stream(msg)
            elif c=='tcl':
                path=self.get_tcldir()
                with open(os.path.join(path,'00PESTIFER-README.txt'),'r') as f:
                    msg=f.read()
                out_stream(msg)

    def get_pdb(self,name):
        return self.pdb_collection.get_pdb(name)

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

    def get_tcldir(self):
        return self.resource_path['tcl']
    
    def get_tcl_pkgdir(self):
        return os.path.join(self.resource_path['tcl'],'pkg')
    
    def get_tcl_scriptsdir(self):
        return os.path.join(self.resource_path['tcl'],'scripts')
    
    def set_pdb_depot(self,depot_path):
        self.pdb_depot_path=depot_path