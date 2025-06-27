# Author: Cameron F. Abrams, <cfa22@drexel.edu>
""" The BaseTask class

"""
import logging
import os
import shutil
import yaml

from .baseobj import BaseObj
from .stringthings import FileCollector

logger=logging.getLogger(__name__)
logging.getLogger("matplotlib").setLevel(logging.WARNING)

class BaseTask(BaseObj):
    """ A base class for Tasks.
    
    Attributes
    ----------
    req_attr: list
        * `specs`: dictionary of specifications; under '`directives`' in yaml
        * `writers`: dictionary of `FileWriters`
        * `prior`: identifier of prior task in sequence
        * `index`: unique integer index of task in run
        * `config`: access to the run config
        * `taskname`: caller-supplied name of task
    
    `subtaskcount`: int
        a count of subtasks in the task
    
    `basename`: string
        a basename used in file-naming conventions for the task

    `statevars`: dict
        dictionary that captures the current state of the run

    `FC`: `FileCollector`
        a structure for tracking files created by the task

    Methods
    -------

    `next_basename(*obj)`:
        determines the next basename based on the task index, subtask
        label (obj[0]) and a subtask count

    `coor_to_pdb()`:
        Generates and executes a VMD script to generate a new PDB
        file from an existing psf and coor file.
    
    `pdb_to_coor()`: 
        Generates and executes a VMD script to generate a new coor
        file from an existing psf and PDB file.
    
    `make_constraint_pdb(dict)`: 
        Generates and executes a VMD script that generates a PDB file 
        with the appropriate attribute-tagging based on the directives 
        in the input dict.

    `inherit_state()`: 
        copies the statevar dict from the previous task onto this
        task's statevar dict

    `save_state(list)`: 
        Based on the current basename, stores the current names of 
        all files whose extension types are in the list
    
    `update_statevars(key,value)`: 
        Updates this task's statevars dict subject to some controls.

    """
    req_attr=BaseObj.req_attr+['specs','config','index','prior','writers','taskname','controller_index']
    yaml_header='generic_task'
    init_msg_options=['INITIATED','STARTED','BEGUN','SET IN MOTION','KICKED OFF','LIT','SPANKED ON THE BOTTOM']

    def __init__(self,config_specs={},controller_specs={}):
        specs=config_specs.copy()
        prior=controller_specs.get('prior',None)
        index=controller_specs.get('index',0)
        writers=controller_specs.get('writers',{})
        taskname=controller_specs.get('taskname','generic_task')
        config=controller_specs.get('config',{})
        controller_index=controller_specs.get('controller_index',0)
        if prior:
            index=prior.index+1
        logger.debug(f'Creating task {taskname} with index {index}')
        input_dict = {
            'index':index,
            'writers':writers,
            'prior':prior,
            'specs':specs,
            'config':config,
            'taskname':taskname,
            'controller_index':controller_index
        }
        super().__init__(input_dict)
        self.subtaskcount=0
        self.statevars={}
        self.FC=FileCollector()
        self.result=0

    def override_taskname(self,taskname):
        logger.debug(f'Overriding task name {self.taskname} to {taskname}')
        self.taskname=taskname

    def do(self):
        return self.result

    def __str__(self):
        return f'{self.controller_index} - {self.index} - {self.taskname} [has prior {self.prior!=None}]'

    def log_message(self,message,**kwargs):
        extra=''
        for k,v in kwargs.items():
            if v:
                extra+=f' ({k}: {v})'
        mtoks=[x.strip() for x in [x.upper() for x in message.split()]]
        if not any([x in self.init_msg_options for x in mtoks]):
            extra+=f' (result: {self.result})'
        logger.info(f'Controller {self.controller_index:02} Task {self.index:02} \'{self.taskname}\' {message} {extra}')

    def get_keepfiles(self):
        """ Returns a list of files that should be kept after the task is done """
        if hasattr(self,'keepfiles'):
            return self.keepfiles
        else:   
            return []
        
    def next_basename(self,*obj):
        label=''
        if len(obj)==1 and len(obj[0])>0:
            label=f'-{obj[0]}'
        default_basename=f'{self.controller_index:02d}-{self.index:02d}-{self.subtaskcount:02d}_{self.taskname}{label}'
        overwrite_basename=self.specs.get('basename',None)
        if overwrite_basename:
            logger.debug(f'Overriding basename {default_basename} with {overwrite_basename}')
            self.basename=overwrite_basename
        else:
            self.basename=default_basename
        self.subtaskcount+=1
    
    def coor_to_pdb(self):
        vm=self.writers['vmd']
        vm.newscript(f'{self.basename}-coor2pdb')
        psf=self.statevars['psf']
        coor=self.statevars['coor']
        vm.addline(f'namdbin2pdb {psf} {coor} {self.basename}.pdb')
        vm.writescript()
        vm.runscript()

    def pdb_to_coor(self):
        vm=self.writers['vmd']
        vm.newscript(f'{self.basename}-pdb2coor')
        pdb=self.statevars['pdb']
        vm.addline(f'pdb2namdbin {pdb} {self.basename}.coor')
        vm.writescript()
        vm.runscript()

    def make_constraint_pdb(self,specs,statekey='consref'):
        vm=self.writers['vmd']
        pdb=self.statevars['pdb']
        force_constant=specs.get('k',self.config['user']['namd']['harmonic']['spring_constant'])
        constrained_atoms_def=specs.get('atoms','all')
        logger.debug(f'constraint spec: {specs["atoms"]}')
        c_pdb=specs.get('consref','')
        if not c_pdb:
            c_pdb=f'{self.basename}-constraints.pdb'
        vm.newscript(f'{self.basename}-make-constraint-pdb')
        vm.addline(f'mol new {pdb}')
        vm.addline(f'set a [atomselect top all]')
        vm.addline(f'$a set occupancy 0.0')
        vm.addline(f'set c [atomselect top "{constrained_atoms_def}"]')
        vm.addline(f'$c set occupancy {force_constant}')
        vm.addline(f'$a writepdb {c_pdb}')
        vm.writescript()
        vm.runscript()
        self.update_statevars(statekey,c_pdb,mode='file')
    
    def inherit_state(self):
        if self.prior:
            logger.debug(f'Inheriting state from prior task {self.prior.taskname}')
            self.statevars=self.prior.statevars.copy()

    def save_state(self,exts=[]):
        for ext in exts:
            mode='strict'
            if ext=='xsc':
                mode='permissive'
            self.update_statevars(ext,f'{self.basename}.{ext}',vtype='file',mode=mode)
        if 'pdb' in exts and 'vel' in self.statevars: # if pdb has changed, assume binary velocity file is stale
            del self.statevars['vel']
        if 'coor' in exts and 'pdb' not in exts:
            self.coor_to_pdb()
            self.update_statevars('pdb',f'{self.basename}.pdb',vtype='file')
        elif 'pdb' in exts and 'coor' not in exts:
            self.pdb_to_coor()
            self.update_statevars('coor',f'{self.basename}.coor',vtype='file')

    def update_statevars(self,key,value,vtype='',mode='strict'):
        if vtype=='file':
            if not os.path.exists(value):
                if mode=='strict':
                    raise FileNotFoundError(f'Expected file "{value}" not found.')
                else:
                    logger.debug(f'Attempt to update state {key} with "{value}" failed since {value} is not a file.')
            else:
                self.statevars[key]=value
        else:
            self.statevars[key]=value

    def write_statefile(self):
        statefilename=self.specs.get('statefile',f'{self.basename}-state.yaml')
        with open(statefilename,'w') as f:
            yaml.dump(self.statevars,f)

    def copy_state(self,exts=[]):
        copies=[]
        for ext in exts:
            if ext in self.statevars:
                ffile=f'{self.basename}.{ext}'
                shutil.copy(self.statevars[ext],ffile)
                self.update_statevars(ext,ffile,vtype='file')
                copies.append(ffile)
        logger.debug(f'Copy state inherited {", ".join(copies)}')