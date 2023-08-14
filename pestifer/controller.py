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

class Task(BaseMod):
    def __init__(self,owner=None,prior=None,next=None):
        input_dict={
            'owner':owner,
            'prior':prior,
            'next':next
        }
        super().__init__(input_dict)

    def do_build(self,bc=None,fc=None):
        pg=Psfgen(self.config.resman,bc)
        pg.beginscript()
        pg.describe_molecule(self.owner.base_molecule,self.mods,file_collector=fc)
        pg.transform_postmods(f'{self.owner.name}',file_collector=fc)
        pg.global_postmods()
        pg.endscript()
        pg.writescript(step.name)
        o,e=pg.runscript()
    
    def do_solvate(self):
        pass

    def do_smdclose(self):
        pass

    def do_relax(self):
        pass


class Step(BaseMod):
    req_attr=BaseMod.req_attr+['name','source','tasks','pdbs']
    def __init__(self,input_dict):
        if not 'pdbs' in input_dict:
            input_dict['pdbs']=[]
        super().__init__(input_dict)
        self.chainIDmanager=ChainIDManager()
        self.tasks_resolved=False
        self.mods_resolved=False

    def resolve_mods(self,mod_classes,modlist_classes):
        if not self.tasks_resolved:
            logger.warning(f'Attempting mod resolution before task resolution')
            return self
        replace_mods={}
        for t in self.tasks:
            if hasattr(t,'build'):
                build_task=t
                mod_list=t.build.get('mods',[])
        if mod_list:
            # any list items that are (a) strings and (b) have '.yaml' or '.yml' at the end
            # are interpreted as files; otherwise, we have a mod_dict
            for m in mod_list:
                if type(m)==str and (m.endswith('.yaml') or m.endswith('.yml')):
                    modfile=m
                    with open(modfile,'r') as f:
                        raw=yaml.safe_load(f)
                        from_file=modparse(raw,mod_classes,modlist_classes)
                        replace_mods=special_update(replace_mods,from_file)
                elif type(m)==dict:
                    mod_dict=m
                    from_userconfig=modparse(mod_dict,mod_classes,modlist_classes)
                    replace_mods=special_update(replace_mods,from_userconfig)
                else:
                    logger.warning(f'Mod list element {m} ignored.')
        build_task['mods']=replace_mods
        for name,mod in build_task['mods'].items():
            if hasattr(mod,'source'):
                if 'pdb' in mod.source:
                    self.pdbs.append(mod.source['pdb'])
        self.mods_resolved=True
        return self
    
    def resolve_tasks(self):
        assert type(self.tasks)==list
        assert len(self.tasks)>1
        assert type(self.tasks[0])==dict
        tasklist=[]
        for tdict in self.tasks:
            tasklist.append(Task(tdict,owner=self))
        self.tasks=tasklist
        self.tasks_resolved=True
        return self

    def do_tasks(self,bc=None,fc=None):
        for task in self.tasks:
            task.do(bc=bc,fc=fc)

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
        self.steps=[]
        self.register_mod_classes()
        self.register_modlist_classes()
        if 'steps' in self.config.defs:
            for step in self.config.defs['steps']:
                self.steps.append(Step(step).resolve_tasks().resolve_mods(self.mod_classes,self.modlist_classes))
    
    def register_mod_classes(self):
        self.mod_classes={}
        for name,cls in inspect.getmembers(sys.modules['pestifer.mods'], lambda x: inspect.isclass(x) and (x.__module__=='pestifer.mods') and 'List' not in x.__name__):
            self.mod_classes[name]=cls

    def register_modlist_classes(self):
        self.modlist_classes={}
        for name,cls in inspect.getmembers(sys.modules['pestifer.mods'], lambda x: inspect.isclass(x) and (x.__module__=='pestifer.mods') and 'List' in x.__name__):
            self.modlist_classes[name]=cls

    def do_steps(self,**kwargs):
        self.check()
        fc=FileCollector()
        bc=ByteCollector()
        for step in self.steps:
            step.injest_molecules()
            step.do_tasks(byte_collector=bc,file_collector=fc)
            if kwargs.get('clean_up',False): fc.flush()

    def check(self):
        stepnames=list(set([x.name for x in self.steps]))
        assert len(stepnames)==len(self.steps),f'Please use unique step names in your BuildSteps section'
