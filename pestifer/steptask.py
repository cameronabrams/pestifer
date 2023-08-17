from .basemod import BaseMod, ModList
from .molecule import Molecule
from .chainids import ChainIDManager
from .util import special_update
import yaml
import os
import logging
logger=logging.getLogger(__name__)
from .util import inspect_classes

class TaskList(ModList):
    def __init__(self,data):
        super().__init__(data)
        self.mod_classes,self.modlist_classes=inspect_classes('pestifer.mods','List')
        self.task_classes=inspect_classes('pestifer.tasks',use_yaml_headers_as_keys=True)
    def append(self,taskdict:dict,owner_step=None):
        taskname=list(taskdict.keys())[0]
        assert taskname in self.task_classes,f'Task "{taskname}" is not recognized; currently implemented tasks are {",".join(self.task_classes)}'
        cls=self.task_classes[taskname]
        task=cls(taskdict)
        task.owner=owner_step
        if len(self)>0:
            task.prior=self[-1]
            task.prior.next=task
        # [x for x in task.__dict__.keys() if not x in task.req_attr][0]
        logging.debug(f'Appending task "{taskname}" to step {task.owner.name}')
        logging.debug(f' -> task specs: {task.specs}')
        if taskname=='terminate_step':
            task.basename=f'{task.owner.name}'
        else:
            task.basename=f'{task.owner.name}-{task.taskname}'
        ''' resolve mods '''
        modlist=task.specs.get('mods',[])
        if modlist:
            for m in modlist:
                if type(m)==str and (m.endswith('.yaml') or m.endswith('.yml')):
                    modfile=m
                    with open(modfile,'r') as f:
                        raw=yaml.safe_load(f)
                        from_file=self.modparse(raw)
                        replace_mods=special_update(replace_mods,from_file)
                elif type(m)==dict:
                    mod_dict=m
                    from_userconfig=self.modparse(mod_dict)
                    replace_mods=special_update(replace_mods,from_userconfig)
                else:
                    logger.warning(f'Mod list element {m} ignored.')
            task.specs['mods']=replace_mods
            for name,mod in task.specs['mods'].items():
                if hasattr(mod,'source'):
                    if 'pdb' in mod.source:
                        self.pdbs.append(mod.source['pdb'])
            task.specs['mods_resolved']=True
        super().append(task)

    def do_tasks(self):
        for task in self:
            task.do()

    def modparse(self,input_dict):
        mod_classes,modlist_classes=self.inspect_classes('pestifer.mods','List')
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
    
class Step(BaseMod):
    req_attr=BaseMod.req_attr+['name','source','tasks','pdbs']
    def __init__(self,input_dict):
        if not 'pdbs' in input_dict:
            input_dict['pdbs']=[]
        super().__init__(input_dict)
        self.chainIDmanager=ChainIDManager()
        self.tasks_resolved=False
    
    def resolve_tasks(self,psfgen,vmdtcl):
        assert type(self.tasks)==list
        assert len(self.tasks)>0
        assert type(self.tasks[0])==dict
        self.psfgen=psfgen
        self.vmdtcl=vmdtcl
        self.tasklist=TaskList([])
        for tdict in self.tasks:
            self.tasklist.append(tdict,owner_step=self)
        self.tasklist.append({'terminate_step':None},owner_step=self)
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

    def do_tasks(self):
        self.tasklist.do_tasks()
