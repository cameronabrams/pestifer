"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import inspect
import sys
import os
import yaml
import logging
logger=logging.getLogger(__name__)

from .config import ConfigSetup
from .molecule import Molecule
from .basemod import BaseMod
from .chainids import ChainIDManager
from .util import special_update
from .psfgen import Psfgen
from .stringthings import FileCollector, ByteCollector

def modparse(input_dict,mod_classes,modlist_classes):
    retdict={}
    for hdr,entries in input_dict.items():
        class_name=[name for name,cls in mod_classes.items() if cls.yaml_header==hdr][0]
        cls=mod_classes[class_name]
        LCls=modlist_classes.get(f'{class_name}List',list)
        for entry in entries:
            assert type(entry) in [dict,str]
            newmod=cls(entry)
            newmod.source='USER'
            if not hdr in retdict:
                retdict[hdr]=LCls([])
            retdict[hdr].append(newmod)
    return retdict

class BuildStep(BaseMod):
    req_attr=BaseMod.req_attr+['name','source','mods','task','pdbs']
    def __init__(self,input_dict):
        if not 'mods' in input_dict:
            input_dict['mods']={}
        if not 'pdbs' in input_dict:
            input_dict['pdbs']=[]
        super().__init__(input_dict)
        self.chainIDmanager=ChainIDManager()

    def resolve_mods(self,mod_classes,modlist_classes):
        replace_mods={}
        mod_dict=self.mods.copy()
        for modfile in mod_dict.get('ModFiles',[]):
            with open(modfile,'r') as f:
                raw=yaml.safe_load(f)
                from_file=modparse(raw,mod_classes,modlist_classes)
                replace_mods=special_update(replace_mods,from_file)
        if 'ModFiles' in mod_dict:
            del mod_dict['ModFiles']
        from_userconfig=modparse(mod_dict,mod_classes,modlist_classes)
        replace_mods=special_update(replace_mods,from_userconfig)
        self.mods=replace_mods
        for name,mod in self.mods.items():
            if hasattr(mod,'source'):
                if 'pdb' in mod.source:
                    self.pdbs.append(mod.source['pdb'])
        return self
    
    def injest_molecules(self):
        self.molecules={}
        if type(self.source)==dict:
            basePDB=self.source.get('pdb',None)
            baseBA=self.source.get('biological_assembly',None)
        elif type(self.source)==str:
            basePDB=self.source
            baseBA=0
        self.molecules[basePDB]=Molecule(source=basePDB).activate_biological_assembly(baseBA,self.chainIDmanager)
        self.base_molecule=self.molecules[basePDB]
        for p in self.pdbs:
            self.molecules[p]=Molecule(source=p)
        return self

class Controller:
    def __init__(self,userconfigfilename):
        self.config=ConfigSetup(userconfigfilename)
        self.build_steps=[]
        self.register_mod_classes()
        self.register_modlist_classes()
        if 'BuildSteps' in self.config.defs:
            for step in self.config.defs['BuildSteps']:
                self.build_steps.append(BuildStep(step).resolve_mods(self.mod_classes,self.modlist_classes))
    
    def register_mod_classes(self):
        self.mod_classes={}
        for name,cls in inspect.getmembers(sys.modules['pestifer.mods'], lambda x: inspect.isclass(x) and (x.__module__=='pestifer.mods') and 'List' not in x.__name__):
            self.mod_classes[name]=cls

    def register_modlist_classes(self):
        self.modlist_classes={}
        for name,cls in inspect.getmembers(sys.modules['pestifer.mods'], lambda x: inspect.isclass(x) and (x.__module__=='pestifer.mods') and 'List' in x.__name__):
            self.modlist_classes[name]=cls

    def do(self,**kwargs):
        self.check()
        fc=FileCollector()
        bc=ByteCollector()
        pg=Psfgen(self.config.resman,bc)
        for step in self.build_steps:
            step.injest_molecules()
            pg.beginscript()
            pg.describe_molecule(step.base_molecule,step.mods,file_collector=fc)
            pg.transform_postmods(step.output,file_collector=fc)
            pg.global_postmods()
            pg.endscript()
            pg.writescript(step.name)
            o,e=pg.runscript()
            if kwargs.get('clean_up',False): fc.flush()

    def check(self):
        stepnames=list(set([x.name for x in self.steps]))
        assert len(stepnames)==len(self.steps),f'Please use unique step names in your BuildSteps section'
