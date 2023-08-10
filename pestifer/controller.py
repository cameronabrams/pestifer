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

def modparse(input_dict,mod_classes):
    retdict={}
    for hdr,entries in input_dict.items():
        Cls=[cls for name,cls in mod_classes.items() if cls.yaml_header==hdr]
        assert len(Cls)==1,f'Cannot resolve {hdr} from available mod classes {[cls.yaml_header for cls in mod_classes.values()]}'
        cls=Cls[0]
        for entry in entries:
            if 'shortcode' in entry:
                newmod=cls.from_shortcode(entry['shortcode'])
            else:
                newmod=cls(entry) # mods by default initialize from dicts
            if not hdr in retdict:
                retdict[hdr]=[]
            retdict[hdr].append(newmod)
    return retdict

class BuildStep(BaseMod):
    req_attr=BaseMod.req_attr+['solvate','mods']
    opt_attr=BaseMod.opt_attr+['smdclose','relax_steps']
    def __init__(self,input_dict):
        if not 'solvate' in input_dict:
            input_dict['solvate']=False
        if not 'mods' in input_dict:
            input_dict['mods']={}
        super().__init__(input_dict)

    def resolve_mods(self,mod_classes):
        replace_mods={}
        mod_dict=self.mods.copy()
        for modfile in mod_dict.get('ModFiles',[]):
            with open(modfile,'r') as f:
                raw=yaml.safe_load(f)
                from_file=modparse(raw,mod_classes)
                replace_mods=special_update(replace_mods,from_file)
        if 'ModFiles' in mod_dict:
            del mod_dict['ModFiles']
        from_userconfig=modparse(mod_dict,mod_classes)
        replace_mods=special_update(replace_mods,from_userconfig)
        self.mods=replace_mods
        return self

class Controller:
    def __init__(self,userconfigfilename):
        self.config=ConfigSetup(userconfigfilename)
        self.chainIDmanager=ChainIDManager()
        self.build_steps=[]
        self.register_mod_classes()
        if 'BuildSteps' in self.config.defs:
            for stage in self.config.defs['BuildSteps']:
                self.build_steps.append(BuildStep(stage).resolve_mods(self.mod_classes))
    
    def register_mod_classes(self):
        self.mod_classes={}
        for name,cls in inspect.getmembers(sys.modules['pestifer.mods'], lambda x: inspect.isclass(x) and (x.__module__=='pestifer.mods') and 'List' not in x.__name__):
            self.mod_classes[name]=cls

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
        return self

    # def report(self):

    #     print(str(self.config))
    #     print(str(self.resman))