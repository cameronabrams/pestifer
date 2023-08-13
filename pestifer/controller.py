"""

.. module:: psfgen
   :synopsis: Manages script generation and execution of psfgen under vmd
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import inspect
import sys
import yaml
import logging
logger=logging.getLogger(__name__)

from .config import ConfigSetup
from .molecule import Molecule
from .basemod import BaseMod
from .chainids import ChainIDManager
from .util import special_update
from .psfgen import Psfgen

def modparse(input_dict,mod_classes,modlist_classes):
    retdict={}
    for hdr,entries in input_dict.items():
        class_name=[name for name,cls in mod_classes.items() if cls.yaml_header==hdr][0]
        cls=mod_classes[class_name]
        LCls=modlist_classes.get(f'{class_name}List',list)
        for entry in entries:
            if 'shortcode' in entry:
                newmod=cls.from_shortcode(entry['shortcode'])
            else:
                entry['source']='USER'
                newmod=cls(entry) # mods by default initialize from dicts
            if not hdr in retdict:
                # print(f'new modlist type {LCls} for class {class_name}')
                # print(f'mod_classes {mod_classes}')
                retdict[hdr]=LCls([])
            retdict[hdr].append(newmod)
    return retdict

class BuildStep(BaseMod):
    req_attr=BaseMod.req_attr+['solvate','mods','output']
    opt_attr=BaseMod.opt_attr+['smdclose','relax_steps']
    def __init__(self,input_dict):
        if not 'solvate' in input_dict:
            input_dict['solvate']=False
        if not 'mods' in input_dict:
            input_dict['mods']={}
        super().__init__(input_dict)

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
        return self
    
    def my_pdbs(self):
        pdbs=[]
        for name,mod in self.mods.items():
            if hasattr(mod,'Source'):
                pdbs.append(mod.Source['pdb'])
        return pdbs

class Controller:
    def __init__(self,userconfigfilename):
        self.config=ConfigSetup(userconfigfilename)
        self.chainIDmanager=ChainIDManager()
        self.build_steps=[]
        self.psfgen=Psfgen(self.config.resman)
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

    def build_molecules(self):
        self.molecules={}
        source_dict=self.config.defs.get('Source',{})
        if not source_dict:
            return self
        basePDB=source_dict.get('pdb','')
        if not basePDB:
            return self
        biological_assembly_index=source_dict.get('biological_assembly',-1)
        self.molecules[basePDB]=Molecule.from_rcsb(pdb_code=basePDB).activate_biological_assembly(biological_assembly_index,self.chainIDmanager)
        # TODO: gather any supplemental molecules
        for step in self.build_steps:
            pdbs_from_step=step.my_pdbs()
            for p in pdbs_from_step:
                if not p in self.molecules:
                    self.molecules[p]=Molecule.from_rcsb(pdb_code=p)
        return self
    
    def do(self):
        self.check()
        self.build_molecules()
        self.psfgen.beginscript()
        for code,mol in self.molecules.items():
            self.psfgen.set_molecule(mol)
        for step in self.build_steps:
            main_pdb=self.config.defs['Source'].get('pdb',None)
            self.psfgen.describe_molecule(self.molecules[main_pdb],step.mods)
            self.psfgen.transform_postmods(step.output)
        self.psfgen.global_postmods()
        self.psfgen.endscript()
        # TODO: run VMD!
    def check(self):
        assert type(self.config.defs)==dict
        assert 'Source' in self.config.defs
        assert 'pdb' in self.config.defs['Source']
        assert self.config.defs['Source'].get('pdb',None)
        assert self.config.defs['Source'].get('biological_assembly',None)

    # def report(self):

    #     print(str(self.config))
    #     print(str(self.resman))