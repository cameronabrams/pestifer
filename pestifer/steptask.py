from .basemod import BaseMod, ModList
from .molecule import Molecule
from .chainids import ChainIDManager
from .util import special_update
import yaml
import os
import sys
import shutil
import inspect
import logging
logger=logging.getLogger(__name__)
    
def inspect_mod_classes():
    mod_classes={}
    for name,cls in inspect.getmembers(sys.modules['pestifer.mods'], lambda x: inspect.isclass(x) and (x.__module__=='pestifer.mods') and 'List' not in x.__name__):
        mod_classes[name]=cls
    modlist_classes={}
    for name,cls in inspect.getmembers(sys.modules['pestifer.mods'], lambda x: inspect.isclass(x) and (x.__module__=='pestifer.mods') and 'List' in x.__name__):
        modlist_classes[name]=cls
    return mod_classes,modlist_classes

def modparse(input_dict):
    mod_classes,modlist_classes=inspect_mod_classes()
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
    def __init__(self,input_dict,owner=None):
        input_dict.update({
            'owner':owner,
            'prior':None,
            'next':None
        })
        assert len(input_dict)==4
        super().__init__(input_dict)

class TaskList(ModList):
    def append(self,task:Task):
        if len(self)>0:
            task.prior=self[-1]
            task.prior.next=task
        assert len(task.__dict__)==4,f'Attempting to append a task that has already been appended?'
        taskname=[x for x in task.__dict__.keys() if not x in ['owner','prior','next']][0]
        assert taskname in ['build','relax','relax_series','smdclose','solvate','terminate_step']
        if   taskname=='build':        taskfunc=self.do_build
        elif taskname=='relax':        taskfunc=self.do_relax
        elif taskname=='relax_series': taskfunc=self.do_relax_series
        elif taskname=='solvate':      taskfunc=self.do_solvate
        elif taskname=='smdclose':     taskfunc=self.do_smdclose
        elif taskname=='terminate_step':     taskfunc=self.terminate_step
        else:
            taskfunc=self.no_func
        task.taskname=taskname
        task.taskfunc=taskfunc
        if taskname=='terminate_step':
            task.basename=f'{self.owner.name}'
        else:
            task.basename=f'{self.owner.name}-{task.taskname}'
        ''' resolve mods '''
        specs=self.__dict__[taskname]        
        modlist=specs.get('mods',[])
        if modlist:
            for m in modlist:
                if type(m)==str and (m.endswith('.yaml') or m.endswith('.yml')):
                    modfile=m
                    with open(modfile,'r') as f:
                        raw=yaml.safe_load(f)
                        from_file=modparse(raw)
                        replace_mods=special_update(replace_mods,from_file)
                elif type(m)==dict:
                    mod_dict=m
                    from_userconfig=modparse(mod_dict)
                    replace_mods=special_update(replace_mods,from_userconfig)
                else:
                    logger.warning(f'Mod list element {m} ignored.')
            specs['mods']=replace_mods
            for name,mod in specs['mods'].items():
                if hasattr(mod,'source'):
                    if 'pdb' in mod.source:
                        self.pdbs.append(mod.source['pdb'])
            specs['mods_resolved']=True
        super().append(task)
        
    def do_build(self,task):
        specs=task.build
        logger.debug(f'do_build: specs {specs}')
        base_molecule=self.owner.base_molecule
        pg=self.psfgen
        pg.beginscript(task.basename)
        pg.topo_aliases()
        pg.describe_molecule(base_molecule,specs['mods'])
        pg.transform_postmods()
        pg.global_postmods()
        pg.endscript()
        pg.writescript()
        o,e=pg.runscript()
        pg.cleanup()
    
    def do_solvate(self,task):
        specs=task.solvate
        logger.debug(f'do_solvate: specs {specs}')
        pass

    def do_smdclose(self,task):
        specs=task.smdclose
        logger.debug(f'do_smdclose: specs {specs}')
        pass

    def do_relax(self,task):
        specs=task.relax
        logger.debug(f'do_relax: specs {specs}')
        pass

    def do_relax_series(self,task):
        specs=task.relax_series
        logger.debug(f'do_relax_series: specs {specs}')
        pass

    def do_terminate_step(self,task):
        prior_basename=task.prior.basename
        logger.debug(f'do_terminate_step')
        for ext in ['.psf','.pdb']:
            shutil.copyfile(f'{prior_basename}{ext}',f'{task.basename}{ext}')

    def do_tasks(self):
        for task in self:
            task.taskfunc(task)

class Step(BaseMod):
    req_attr=BaseMod.req_attr+['name','source','tasks','pdbs']
    def __init__(self,input_dict):
        if not 'pdbs' in input_dict:
            input_dict['pdbs']=[]
        super().__init__(input_dict)
        self.chainIDmanager=ChainIDManager()
        self.tasks_resolved=False
    
    def resolve_tasks(self):
        assert type(self.tasks)==list
        assert len(self.tasks)>1
        assert type(self.tasks[0])==dict
        self.tasklist=TaskList([])
        for tdict in self.tasks:
            self.tasklist.append(Task(tdict,owner=self))
        self.tasklist.append(Task({'terminate_step':None},owner=self))
        self.tasks_resolved=True
        return self

    def injest_molecules(self):
        self.molecules={}
        psf_exists=False
        if type(self.source)==dict:
            basePDB=self.source.get('pdb',None)
            baseBA=self.source.get('biological_assembly',None)
        elif type(self.source)==str:
            basePDB=self.source
            assert os.path.exists(f'{basePDB}.pdb')
            psf_exists=os.path.exists(f'{basePDB}.psf')
            baseBA=0
        self.molecules[basePDB]=Molecule(source=basePDB,use_psf=psf_exists).activate_biological_assembly(baseBA,self.chainIDmanager)
        self.base_molecule=self.molecules[basePDB]
        for p in self.pdbs:
            self.molecules[p]=Molecule(source=p)
        return self

